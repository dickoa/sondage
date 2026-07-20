#ifndef SONDAGE_SAMPLING_CORE_H
#define SONDAGE_SAMPLING_CORE_H

/* Count all probabilities above threshold, but never write beyond capacity. */
static int sampling_extract_selected(const double *prob, int N,
                                     double threshold, int *out,
                                     int capacity) {
    int selected = 0;
    for (int i = 0; i < N; i++) {
        if (prob[i] > threshold) {
            if (selected < capacity) out[selected] = i + 1;
            selected++;
        }
    }
    return selected;
}

#endif
