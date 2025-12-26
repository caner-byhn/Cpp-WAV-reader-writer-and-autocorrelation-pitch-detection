#include <iostream>
#include <cstdint>
#include <vector>
#include <fstream>
#include <limits>
#include <algorithm>
#include <math.h>


constexpr float INT16_MAX_F = 
    static_cast<float>(std::numeric_limits<int16_t>::max()); // this is 32767, int16 has the limits of -32768 to 32767

constexpr float INT16_SCALE =
    -static_cast<float>(std::numeric_limits<int16_t>::min()); // this is 32768


struct WavHeader {
    char riff[4] = {'R','I','F','F'};
    uint32_t chunkSize;
    char wave[4] = {'W','A','V','E'};

    char fmt[4] = {'f','m','t',' '};
    uint32_t fmtSize = 16;
    uint16_t audioFormat = 1;
    uint16_t numChannels;
    uint32_t sampleRate;
    uint32_t byteRate;
    uint16_t blockAlign;
    uint16_t bitsPerSample = 16;

    char data[4] = {'d','a','t','a'};
    uint32_t dataSize;
};


void writeWav(
    const std::string& filename,
    const std::vector<float>& samples,
    uint32_t sampleRate,
    uint16_t numChannles);

int16_t floatToInt16(float s);
std::vector<float> mixToMono(const std::vector<float>& samples, uint16_t channels);
float autocorrelation(const std::vector<float>& samples, uint32_t sampleRate);


int main() {
    std::vector<float> samples;
    WavHeader header;
    size_t count = 0;

    std::ifstream file("D.wav", std::ios::binary);
    file.read(reinterpret_cast<char*>(&header), sizeof(header));

    size_t totalSamples = header.dataSize / (header.bitsPerSample / 8);
    samples.reserve(totalSamples);

    int16_t sample;
    while (file.read(reinterpret_cast<char*>(&sample), sizeof(sample))) {
        samples.push_back(sample / INT16_SCALE);
        if (++count % 1'000'000 == 0) {
            std::cout << count << " samples\n";
        }
    }

    for (size_t i = 0; i < samples.size(); i += 2) {
        float step = 1.0f / 256.0f;
        float dither = (rand() / float(RAND_MAX) - 0.5f) * step;

        samples[i]     = std::round((samples[i] + dither) / step) * step;
        samples[i + 1] = std::round((samples[i + 1] + dither) / step) * step;
    }

    std::vector<float> samplesLeft;
    samplesLeft.reserve(totalSamples / header.numChannels);

    for (size_t i = 0; i < samples.size(); i+= header.numChannels) {
        samplesLeft.push_back(samples[i]);
    }

    std::vector<float> mono = mixToMono(samples, header.numChannels);

    float hzStereo = autocorrelation(samplesLeft, header.sampleRate);
    float hzMono = autocorrelation(mono, header.sampleRate);

    std::cout << "Channels: " << header.numChannels << "\n";
    std::cout << "Stereo single channel: " << hzStereo << "\n";
    std::cout << "Mono converted: " << hzMono << "\n";
    
    float note = 12.0f * log2(hzMono / 440.0f) + 69.0f;
    std::cout << "Note: " << note << "\n";
    int nearest = std::round(note);
    float cents = (note - nearest) * 100.0f;

    /*
    if (cents > 5)
        "Tune down"
    else if (cents < -5)
        "Tune up"
    else
        "In tune"
    */

    writeWav(
        "output.wav",
        samples,
        header.sampleRate,
        header.numChannels);
}


void writeWav(
    const std::string& filename,
    const std::vector<float>& samples,
    uint32_t sampleRate,
    uint16_t numChannles
) {
    WavHeader header;

    header.numChannels = numChannles;
    header.sampleRate = sampleRate;
    header.bitsPerSample = 16;
    header.blockAlign = numChannles * sizeof(int16_t);
    header.byteRate = sampleRate * header.blockAlign;
    header.dataSize = samples.size() * sizeof(int16_t);
    header.chunkSize = 36 + header.dataSize;

    std::ofstream out(filename, std::ios::binary | std::ios::out | std::ios::trunc);

    out.write(reinterpret_cast<char*>(&header), sizeof(header));

    for (float s : samples) {
        int16_t pcm = floatToInt16(s);
        out.write(reinterpret_cast<char*>(&pcm), sizeof(pcm));
    }
}


int16_t floatToInt16(float s) {
    s = std::clamp(s, -1.0f, 1.0f);
    return static_cast<int16_t>( std::round(s * INT16_MAX_F) );
}


std::vector<float> mixToMono(const std::vector<float>& samples, uint16_t channels) {
    std::vector<float> mono;
    size_t frames = samples.size() / channels;
    mono.reserve(frames);

    if (samples.size() % channels != 0) {
        std::cout << "Invalid audio file\n";
        return mono;
    }

    for (size_t i = 0; i < frames; i += channels) {
        float sum = 0.0f;

        for (int c = 0; c < channels; c++) {
            sum += samples[i + c];
        }
        mono.push_back(sum / channels);
    }

    return mono;
}


float autocorrelation(const std::vector<float>& samples, uint32_t sampleRate) {
    float minFreq = 60.0f;
    float maxFreq = 500.0f;

    int minLag = static_cast<int>(sampleRate / maxFreq);
    int maxLag = static_cast<int>(sampleRate / minFreq);

    float bestLag = 0.0f;
    float bestCorr = -1.0f;

    for (int lag = minLag; lag <= maxLag; lag++) {

        float dot = 0.0f;
        float energy1 = 0.0f;
        float energy2 = 0.0f;

        for (int i = 0; i + lag < samples.size(); i++) {

            float a = samples[i];
            float b = samples[i + lag];

            dot += a * b;
            energy1 += a * a;
            energy2 += b * b;
        }
        
        if (energy1 == 0.0f || energy2 == 0.0f) continue;

        float normalized = dot / (sqrt(energy1) * sqrt(energy2));
        if (normalized > bestCorr) {
            bestCorr = normalized;
            bestLag = lag;
        }
    }

    return bestLag > 0 ? static_cast<float>(sampleRate) / bestLag : 0.0f;
}