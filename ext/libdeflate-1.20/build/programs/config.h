#ifndef CONFIG_H
#define CONFIG_H

/* Is the clock_gettime() function available? */
#define HAVE_CLOCK_GETTIME

/* Is the futimens() function available? */
#define HAVE_FUTIMENS

/* Is the posix_fadvise() function available? */
#define HAVE_POSIX_FADVISE

/* Is the posix_madvise() function available? */
#define HAVE_POSIX_MADVISE

/* Does stat() provide nanosecond-precision timestamps? */
#define HAVE_STAT_NANOSECOND_PRECISION

#endif /* CONFIG_H */
