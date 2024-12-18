#ifndef ATK_TIMING_H
#define ATK_TIMING_H

#include <time.h>
#include <sys/time.h>

namespace at
{
    static inline double diffclock(clock_t clock1,clock_t clock2)
    {
      double diffticks=clock2-clock1;
      return diffticks/CLOCKS_PER_SEC;
    }

    static inline double diffclock(clock_t clock1)
    {
      clock_t clock2 = clock();
      return diffclock(clock1,clock2);
    }
};

#endif /* ATK_TIMING_H */

