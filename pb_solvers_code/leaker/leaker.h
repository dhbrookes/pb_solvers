/*
 * Leaker.h - defines for the Leaker memory leak detection package
 *
 * Copyright 2011 by Dara Hazeghi
 *
 * See README.txt or http://http://left404.com/leaker/ for more information.
 * Special thanks to Marshall Thomas for many of the ideas contained therein.
 *
 * This program is distributed under the terms of the GNU GPL version 2.
 */
/* 0.1 2011-10-15 */

#ifndef _LEAKER_H
#define _LEAKER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef __cplusplus
#include <memory>
#endif

/* Visual C++ fixes */
#ifdef _MSC_VER
/* __func__ is not defined in VC++ */
#define __func__ __FUNCTION__
/* disable warnings given for using standard C functions */
#pragma warning (disable : 4996)
#endif

#define START_SIZE  128     /* Initial size of memory allocation table */

#define GUARD_SIZE  4       /* Padding at the end of each allocated block */
#define GUARD_STR   "\014\033\014\033"  /* magic string to pad allocation */

typedef struct _LEAK_T
{
    void *addr;             /* address of memory allocated              */
    const char *alloc;      /* type of allocator used                   */
    size_t size;            /* size of allocated memory (bytes)         */
    size_t sequence;		/* relative position of allocation in code  */
    const char *file;       /* name of file where allocation made       */
    const char *func;       /* name of function where allocation made   */
    size_t line;			/* line number where allocation made        */
    struct _LEAK_T *next;   /* linked list pointer                      */
} _LEAK_T;

typedef struct
{
    _LEAK_T **table;
    size_t rows;			/* number of rows in table                */
    size_t count;			/* number of entries in table             */
    size_t bytes;			/* number of bytes currently allocated    */
    size_t serial;			/* number of next insertion               */
    
    size_t overflows;		/* number of incorrect deallocations      */
    size_t mismatches;		/* number of mismatched allocs/deallocs   */
    size_t bad_frees;		/* number of bad attempts to free         */
} _HTABLE_T;

extern _HTABLE_T _leaker;

/* report information on current memory allocations */
void _Leaker_Dump(void);

/* replacement for standard C allocation and deallocation functions */
void *_malloc(size_t size, const char *file, const char *func,
                     unsigned long line);
void *_calloc(unsigned long count, size_t size, const char *file,
                     const char *func, unsigned long line);
void *_realloc(void *ptr, size_t size, const char *file,
                      const char *func, unsigned long line);
void _free(void *ptr, const char *file, const char *func,
                     unsigned long line);

/* preprocessor magic to override built-in allocation functions with our own */
#define malloc(size)		_malloc(size, __FILE__, __func__, __LINE__)
#define calloc(n, size)		_calloc(n, size, __FILE__, __func__, __LINE__)
#define free(ptr)			_free(ptr, __FILE__, __func__, __LINE__)
#define realloc(ptr, size)	_realloc(ptr, size, __func__, __FILE__, __LINE__)

#ifdef __cplusplus

/* hackish solution to the problem of overriding C++ operator new/delete */
extern const char *_leaker_file;
extern const char *_leaker_func;
extern unsigned long _leaker_line;

#define new (_leaker_file=__FILE__, _leaker_func=__func__, \
    _leaker_line=__LINE__) && 0 ? NULL : new
#define delete _leaker_file=__FILE__, _leaker_func=__func__, \
    _leaker_line=__LINE__, delete

#endif

#endif
