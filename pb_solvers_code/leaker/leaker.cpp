/*
 * Leaker.c - the Leaker memory leak detection package
 *
 * Copyright 2011 by Dara Hazeghi
 *
 * See README.txt or http://http://left404.com/leaker/ for more information.
 
 * Special thanks to Marshall Thomas for many of the ideas implemented here.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 */
/* 0.1 2011-10-15 */

#ifdef __LEAK
#include "leaker.h"
#endif

/* global table containing allocation information */
_HTABLE_T _leaker = { NULL, 0, 0, 0, 0, 0, 0, 0 };


/* disable macros for internal use */
#undef malloc
#undef calloc
#undef realloc
#undef free

#ifdef __cplusplus
const char *_leaker_file = "unknown";
const char *_leaker_func = "unknown";
unsigned long _leaker_line = 0;

#undef new
#undef delete
#endif

/* internal function prototypes */
static void _Leaker_Init(void);
static void _Leaker_Grow(void);

static void _Leaker_Add(void *addr, size_t size, const char *alloc,
                        const char *file, const char *func, size_t line);
static size_t _Leaker_Remove(void *addr, const char *dealloc, const char *file,
                            const char *func, size_t line);
static _LEAK_T **_Leaker_Find(void *addr);

static unsigned long _Leaker_Hash(void *addr);

static int _Leaker_Check_Dealloc(const char *alloc, const char *dealloc);
static void _Leaker_Init_Guard(void *addr, size_t size);
static int _Leaker_Check_Guard(void *addr, size_t size);
static void _Leaker_Scribble(void *ptr, size_t size);

static void _Leaker_Print_Entry(_LEAK_T *entry);
static void _Leaker_Report(void);

static void _Leaker_Dump_Entry(_LEAK_T *entry);
static _LEAK_T **_Leaker_Build_List(void);

static int _Leaker_Compare_Entries(const void *first, const void *second);

/* dump current allocation information and statistics */
void _Leaker_Dump(void)
{
    size_t overflows = _leaker.overflows;
    
	fprintf(stderr, "\nLeaker report:\n");
    
	if(_leaker.count == 0)
    {
		fprintf(stderr, "No allocations.\n\n");
        return;
    }
    else
    {
        _LEAK_T **table = _Leaker_Build_List();
        unsigned int i;
        fprintf(stderr, "%lu allocations (%lu bytes) in table of %lu rows.\n",
                _leaker.count, _leaker.bytes - _leaker.count * GUARD_SIZE,
                _leaker.rows);
        
        for(i = 0; i < _leaker.count; i++)
        {
            _Leaker_Dump_Entry(table[i]);
        }
        fprintf(stderr, "\n");
        
        free(table);
    }
    
    if(_leaker.mismatches)
        fprintf(stderr, "Mismatches: %lu allocation/deallocations don't match!\n",
                _leaker.mismatches);
    if(_leaker.overflows)
        fprintf(stderr, "Overflows: %lu allocations overflowed (wrote off end)!\n",
                _leaker.overflows);
    if(_leaker.bad_frees)
        fprintf(stderr, "Bad deallocs: %lu attempts made to deallocate unallocated pointers!\n",
                _leaker.bad_frees);
    
    _leaker.overflows = overflows; /* restore current overflow count */
}

/* replacement for malloc */
void *_malloc(size_t size, const char *file, const char *func,
                     unsigned long line)
{
    void *ptr;
    size += GUARD_SIZE;
    
    if(!(ptr = malloc(size)))
    {
        fprintf(stderr,"%s:%s():%i aborting: malloc() failed!\n", __FILE__,
                __func__, __LINE__);
        exit(2);
    }
    
    _Leaker_Init_Guard(ptr, size);
    _Leaker_Add(ptr, size, "malloc", file, func, line);
    return ptr;
}

/* replacement for calloc */
void *_calloc(unsigned long count, size_t size, const char *file,
                     const char *func, unsigned long line)
{
    void *ptr;
    size = size * count + GUARD_SIZE;
    
    if(!(ptr = calloc(1, size)))
    {
        fprintf(stderr,"%s:%s():%i aborting: calloc() failed!\n", __FILE__,
                __func__, __LINE__);
        exit(2);
    }
    
    _Leaker_Init_Guard(ptr, size);
    _Leaker_Add(ptr, size, "calloc", file, func, line);
    return ptr;
}

/* replacement for realloc */
void *_realloc(void *ptr, size_t size, const char *file,
                      const char *func, unsigned long line)
{
    void *ptr_new;
    size_t old_size = 0, len;
    
    size += GUARD_SIZE;
    
    /* check if pointer given to realloc is valid */
    if(ptr && !(old_size = _Leaker_Remove(ptr, "realloc", file, func, line)))
    {
        ptr = NULL;
    }
    
    /* to help catch realloc errors, ensure each realloc is at a new address */
    if(!(ptr_new = malloc(size)))
    {
        fprintf(stderr,"%s:%s():%i aborting: realloc() failed!\n",
                __FILE__, __func__, __LINE__);
        exit(2);
    }
    
    /* if realloc was given a valid pointer, copy over the old data */
    if(ptr)
    {
        if(old_size > size) len = size - GUARD_SIZE;
        else len = old_size - GUARD_SIZE;
    
        memcpy(ptr_new, ptr, len);
        
        _Leaker_Scribble(ptr, old_size);
        free(ptr);
    }
        
    _Leaker_Init_Guard(ptr_new, size);
    _Leaker_Add(ptr_new, size, "realloc", file, func, line);
    return ptr_new;
}

/* replacement for free */
void _free(void *ptr, const char *file, const char *func,
                     unsigned long line)
{
    size_t size;
    if((size = _Leaker_Remove(ptr, "free", file, func, line)))
    {
        _Leaker_Scribble(ptr, size);
        free(ptr);
    }
}

#ifdef __cplusplus

/* replacement for operator new */
void* operator new (size_t size) throw (std::bad_alloc)
{
    void *ptr = NULL;
    size += GUARD_SIZE;
    
    if(!(ptr = malloc(size)))
    {
        fprintf(stderr,"%s:%s():%i aborting: calloc() failed!\n",
                __FILE__, __func__, __LINE__);
        exit(2);
    }
    
    _Leaker_Init_Guard(ptr, size);
    _Leaker_Add(ptr, size, "new", _leaker_file, _leaker_func, _leaker_line);

    /* in case new is called from library code where macro has not overriden
     * new and updated _leaker_file, _leaker_func, etc */
    _leaker_file = "unknown";
    _leaker_func = "unknown";
    _leaker_line = 0;

    return ptr;
}

/* replacement for operator vector new */
void* operator new [] (size_t size) throw (std::bad_alloc)
{
    void *ptr;
    size += GUARD_SIZE;
    
    if(!(ptr = malloc(size)))
    {
        fprintf(stderr,"%s:%s():%i aborting: calloc failed!\n",
                __FILE__, __func__, __LINE__);
        exit(2);
    }
    
    _Leaker_Init_Guard(ptr, size);
    _Leaker_Add(ptr, size, "new[]", _leaker_file, _leaker_func, _leaker_line);
    
    _leaker_file = "unknown";
    _leaker_func = "unknown";
    _leaker_line = 0;
    
    return ptr;
}

/* replacement for operator delete */
void operator delete(void* ptr) throw ()
{
    size_t size;
    if((size = _Leaker_Remove(ptr, "delete", _leaker_file, _leaker_func,
                              _leaker_line)))
    {
        _Leaker_Scribble(ptr, size);
        free(ptr);
    }
}

/* replacement for operator vector delete */
void operator delete [] (void* ptr) throw ()
{
    size_t size;
    if((size = _Leaker_Remove(ptr, "delete[]", _leaker_file, _leaker_func,
                              _leaker_line)))
    {
        _Leaker_Scribble(ptr, size);
        free(ptr);
    }    
}

#endif

/* end of user-visible functions */

/* initialize table */
static void _Leaker_Init(void)
{
    _leaker.table = (_LEAK_T **) calloc(START_SIZE, sizeof(_LEAK_T *));
    
    if(!_leaker.table)
    {
        fprintf(stderr,"%s:%s():%i aborting: calloc() for table failed!\n",
                __FILE__, __func__, __LINE__);
        exit(2);
    }
    
    _leaker.rows = START_SIZE;
    
    /* register so that leak information always displayed upon termination */
    atexit(_Leaker_Report);
}

/* grow table to accommodate more entries (by factor of 4) */
static void _Leaker_Grow(void)
{
    unsigned int i;
    
    _HTABLE_T old = _leaker;
    _leaker.rows *= 4;
    _leaker.count = 0;
    _leaker.bytes = 0;
    
    _leaker.table = (_LEAK_T **) calloc(_leaker.rows, sizeof(_LEAK_T *));
    
    if(!_leaker.table)
    {
        fprintf(stderr,"%s:%s():%i aborting: calloc() for table failed!\n",
                __FILE__, __func__, __LINE__);
        exit(2);
    }    
    
    /* remove every entry from the old table and add it to the new */
    for(i = 0; i < old.rows; i++)
    {
        _LEAK_T *mover = old.table[i];
        while(mover)
        {
            _LEAK_T *temp = mover;
            mover = mover->next;
            _leaker.serial = temp->sequence;
            _Leaker_Add(temp->addr, temp->size, temp->alloc,
                        temp->file, temp->func, temp->line);
            free(temp);
        }
    }
    _leaker.serial = old.serial;
    
    free(old.table);
}

/* add a new allocation to the table */
void _Leaker_Add(void *addr, size_t size, const char *alloc,
                        const char *file, const char *func, size_t line)
{
    _LEAK_T **mover, *temp;
    
    if(!_leaker.table) _Leaker_Init();
    if(_leaker.count > _leaker.rows / 2) _Leaker_Grow();
    
    mover = &(_leaker.table[_Leaker_Hash(addr)]);
    while(*mover)
    {
        if((*mover)->addr == addr) break;
        
        mover = &((*mover)->next);
    }
    
    if(*mover) /* if an address is allocated twice, we are in trouble! */
    {
        fprintf(stderr, "%s:%s():%lu fatal error: address %p already in use!\n",
                file, func, line, addr);
        exit(2);
    }
    
    if(!(temp = (_LEAK_T *)malloc(sizeof(_LEAK_T))))
    {
        fprintf(stderr,"%s:%s():%i aborting: malloc() for entry failed!\n",
                __FILE__, __func__, __LINE__);
        exit(2);
    }
    
    temp->addr = addr;
    temp->size = size;
    temp->sequence = _leaker.serial++;
    temp->alloc = alloc;
    temp->file = file;
    temp->func = func;
    temp->line = line;
    temp->next = NULL;
    
    *mover = temp;
    
    _leaker.count++;
    _leaker.bytes += size;
}

/* remove an allocation from the table, and report any inconsistencies */
size_t _Leaker_Remove(void *addr, const char *dealloc, const char *file,
                            const char *func, size_t line)
{
    _LEAK_T **mover, *temp;
    size_t size;
    
    if(!_leaker.table) _Leaker_Init();
    
    mover = _Leaker_Find(addr);
    
    if(!*mover) /* given a bad pointer */
    {
        fprintf(stderr,"\nLEAKER: %s:%s():%lu %s error: pointer was not allocated!\n\n",
                file, func, line, dealloc);
        _leaker.bad_frees++;
        return 0;
    }
    
    temp = *mover;
    *mover = temp->next;
    _leaker.count--;
    _leaker.bytes -= temp->size;
    
    size = temp->size;
        
    if(!_Leaker_Check_Guard(addr, temp->size)) /* guard overwritten */
    {
        fprintf(stderr, "\nLEAKER: %s:%s():%lu checking error: wrote off end of memory allocated at %s:%s():%lu.\n\n",
                file, func, line, temp->file, temp->func, temp->line);
        _leaker.overflows++;
    }
    if(!_Leaker_Check_Dealloc(temp->alloc, dealloc)) /* wrong dealloc function */
    {
        fprintf(stderr, "\nLEAKER: %s:%s():%lu mismatch error: memory allocated at %s:%s():%lu with %s, deallocated with %s.\n\n",
                file, func, line, temp->file, temp->func, temp->line,
                temp->alloc, dealloc);
        _leaker.mismatches++;
    }
        
    free(temp);
    return size;
}

/* find a pointer in the table, return a pointer to a pointer to it */
static _LEAK_T **_Leaker_Find(void *addr)
{
    _LEAK_T **mover = &(_leaker.table[_Leaker_Hash(addr)]);
    while(*mover)
    {
        if((*mover)->addr == addr) break;
        
        mover = &((*mover)->next);
    }
    return mover;
}

/* Thomas Wang's 64-bit hash function - works well for integers, and is
 * significantly faster than the DJB function since.  It is also slightly
 * better in distributing keys.
 * http://www.concentric.net/~Ttwang/tech/inthash.htm
 */
static unsigned long _Leaker_Hash(void *addr)
{
    unsigned long address = (unsigned long) addr;
	address = (~address) + (address << 21); /* (a << 21) - a - 1; */
	address = address ^ (address >> 24);
	address = (address + (address << 3)) + (address << 8); /* a * 265 */
	address = address ^ (address >> 14);
	address = (address + (address << 2)) + (address << 4); /* a * 21 */
	address = address ^ (address >> 28);
	address = address + (address << 31);
	return address % _leaker.rows;
}

/* return 1 if the allocator and deallocator are compatible, 0 otherwise */
static int _Leaker_Check_Dealloc(const char *alloc, const char *dealloc)
{
    if((strcmp(dealloc, "free") == 0) || (strcmp(dealloc, "realloc") == 0))
    {
        if     (strcmp(alloc,  "malloc") == 0) return 1;
        else if(strcmp(alloc,  "calloc") == 0) return 1;
        else if(strcmp(alloc, "realloc") == 0) return 1;
        else return 0;
    }
    if((strcmp(alloc, "new"  ) == 0) && (strcmp(dealloc, "delete"  ) == 0))
        return 1;
    if((strcmp(alloc, "new[]") == 0) && (strcmp(dealloc, "delete[]") == 0))
        return 1;
    
    return 0;
}

/* initialize the guard at the end of the allocation */
static void _Leaker_Init_Guard(void *addr, size_t size)
{
    strncpy((char *)addr + size - GUARD_SIZE, GUARD_STR, GUARD_SIZE);
}

/* verify the guard at the end of the allocation, return 1 if successful */
static int _Leaker_Check_Guard(void *addr, size_t size)
{    
    return (strncmp((char *)addr + size - GUARD_SIZE, GUARD_STR, GUARD_SIZE)
                == 0);
}

/* zap a give allocation - mainly useful to make later mistakes easier to see */
static void _Leaker_Scribble(void *ptr, size_t size)
{
    memset(ptr, '\0', size);
}

/* report leaks and errors, and deallocate all remaining memory */
static void _Leaker_Report(void)
{
    if(!(_leaker.count || _leaker.mismatches || _leaker.overflows
         || _leaker.bad_frees))
    {
        if(_leaker.table) free(_leaker.table);
        return;
    }
    
    fprintf(stderr, "\nLEAKER: errors found!\n");
    
    if(_leaker.count) /* print out list of leaks, clean up */
    {
        unsigned int i;
        _LEAK_T **table = _Leaker_Build_List();
        
        fprintf(stderr, "Leaks found: %lu allocations (%lu bytes).\n",
                _leaker.count, _leaker.bytes - _leaker.count * GUARD_SIZE);
        
        for(i = 0; i < _leaker.count; i++)
        {
            _Leaker_Print_Entry(table[i]);
            free(table[i]->addr);
            free(table[i]);
        }
        
        fprintf(stderr, "\n");
        free(table);
    }
    
	if(_leaker.table) free(_leaker.table);
    
    /* report other errors */
    if(_leaker.mismatches)
        fprintf(stderr, "Mismatches: %lu allocation/deallocations don't match.\n",
                _leaker.mismatches);
    
    if(_leaker.overflows)
        fprintf(stderr, "Overflows: %lu allocations overflowed (wrote off end).\n",
                _leaker.overflows);
    if(_leaker.bad_frees)
        fprintf(stderr, "Bad deallocs: %lu attempts made to deallocate unallocated pointers.\n",
                _leaker.bad_frees);
    
}

/* print the given entry */
static void _Leaker_Print_Entry(_LEAK_T *entry)
{
    fprintf(stderr, "%s:%s():%lu memory leak: memory was not deallocated.\n",
            entry->file, entry->func, entry->line);
    if(!_Leaker_Check_Guard(entry->addr, entry->size))
    {
        fprintf(stderr, "%s:%s():%lu checking error: wrote off end of allocation.\n",
                entry->file, entry->func, entry->line);
        _leaker.overflows++;
    }
}

/* dump the given entry */
static void _Leaker_Dump_Entry(_LEAK_T *entry)
{
    fprintf(stderr, "%s:%s():%lu address: %p bytes: %lu",
            entry->file, entry->func, entry->line, entry->addr,
            entry->size - GUARD_SIZE);
    if(!_Leaker_Check_Guard(entry->addr, entry->size))
    {
        fprintf(stderr, " OVERFLOWED.\n");
        _leaker.overflows++;
    }
    else
        fprintf(stderr, ".\n");
}

/* compare two leak entries based upon their sequence number */
static int _Leaker_Compare_Entries(const void *first, const void *second)
{
    const _LEAK_T *f = *(_LEAK_T **)first, *s = *(_LEAK_T **)second;
    
    return (int)(f->sequence - s->sequence);
}

/* build and return a sorted list of leak entry pointers */
static _LEAK_T **_Leaker_Build_List(void)
{
    _LEAK_T **table, **mover;
    unsigned int i;
    
    if(!_leaker.table || !_leaker.count) return NULL;
    
    if(!(table = (_LEAK_T **)malloc(sizeof(_LEAK_T *) * _leaker.count)))
    {
        fprintf(stderr,"%s:%s():%i aborting: malloc failed!\n",
                __FILE__, __func__, __LINE__);
        exit(2);
    }
    
    mover = table;
    
    for(i = 0; i < _leaker.rows; i++)
    {
        _LEAK_T *temp = _leaker.table[i];
        
        while(temp)
        {
            *mover = temp;
            mover++;
            temp = temp->next;
        }
    }
    
    qsort((void *)table, _leaker.count, sizeof(_LEAK_T *), _Leaker_Compare_Entries);
    
    return table;
}

