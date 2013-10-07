// Copyright 2013, Andrew Ross
// Distributable under the GNU LGPL v2.1 , see COPYING for details
#ifndef _THREAD_HPP
#define _THREAD_HPP

#include <pthread.h>

// Trivial wrapper around basic synchronization primitives in the
// expectation of porting to something other than POSIX.

struct mutex {
    mutex() { pthread_mutex_init(&m, 0); }
    ~mutex() { pthread_mutex_destroy(&m); }
    void lock() { pthread_mutex_lock(&m); }
    void unlock() { pthread_mutex_unlock(&m); }
protected:
    pthread_mutex_t m;
};

struct cvar : public mutex {
    cvar() { pthread_cond_init(&c, 0); }
    ~cvar() { pthread_cond_destroy(&c); }
    void wait() { pthread_cond_wait(&c, &m); }
    void signal() { pthread_cond_signal(&c); }
private:
    pthread_cond_t c;
};

struct sem {
    sem() : n(0) {}
    void up() { cv.lock(); n++; cv.signal(); cv.unlock(); }
    void down() { cv.lock(); while(!n) cv.wait(); n--; cv.unlock(); }
private:
    int n;
    cvar cv;
};

typedef void (*thread_fn)(void*);

inline void start_thread(thread_fn fn, void* arg)
{
    pthread_t t;
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    pthread_create(&t, &attr, (void*(*)(void*))fn, arg);
}

#endif // _THREAD_HPP
