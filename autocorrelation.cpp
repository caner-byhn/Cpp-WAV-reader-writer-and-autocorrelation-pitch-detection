#include <iostream>
#include <cstdint>
#include <vector>
#include <math.h>


// Take a window of samples.
// For each lag > 0, compare the window with a time-shifted version of itself.
// For each lag, compute the sum of sample[i] * sample[i + lag].
// Each sum is a similarity score for that lag.
// All scores together form the autocorrelation curve.
// Peaks in this curve indicate repeating periods in the signal.

// We want to normalize the sum


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


//Use C style ring buffer for processing real time samples.
//Ring buffer is faster than C++ dynamic containers and more compatible with embedded hardware.
/*
constexpr int N = 9;
float history[N];
int index = 0;
int count = 0;

void push(float hz) {
    history[index] = hz;
    index = (index + 1) % N;
    count = std::min(count + 1, N);
}
*/