#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
main()
{
     struct timeval tv;
     /*struct timezone tz;*/
     struct tm *tm;
     gettimeofday(&tv, NULL);
     long int itime;
     itime = (long) tv.tv_usec;
     /*tm=localtime(&tv.tv_sec);*/
     printf(" %ld \n", itime);
     /*printf(" %d:%02d:%02d %d \n", tv->tm_hour, tv->tm_min,
       tm->tm_sec, tv.tv_usec);*/
     exit(0);
}
