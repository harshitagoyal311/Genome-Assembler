#include "data.h"

typedef data* (*samplingFunc)(dataArray*, int);

data *reservoirSampling2(dataArray *myArrayPointer, int k);