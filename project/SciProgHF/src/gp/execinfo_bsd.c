#include "trace.h"

#if defined (__BSD_EXECINFO__)
/* 
   if this is not a linux system, we use the implementation below (from freebsd). 
   note: the copyright notice her should be kept intact!
*/

/*
 * Copyright (c) 2003 Maxim Sobolev <sobomax@FreeBSD.org>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 *  Id: execinfo.c,v 1.3 2004/07/19 05:21:09 sobomax Exp  
 */

#include <sys/types.h>
#include <sys/uio.h>
#include <dlfcn.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "execinfo_bsd.h"
#include "stacktraverse_bsd.h"

#define D10(x) ceil(log10(((x) == 0) ? 2 : ((x) + 1)))

inline static void *
realloc_safe(void *ptr, size_t size)
{
    void *nptr;

    nptr = realloc(ptr, size);
    if (nptr == NULL)
        free(ptr);
       
    return nptr;
}

size_t
backtrace(void **buffer, size_t size)
{
    size_t i;

    for (i = 1; getframeaddr(i + 1) != NULL && i != size + 1; i++) {
        buffer[i - 1] = getreturnaddr(i);
        if (buffer[i - 1] == NULL)
            break;
    }

    return i - 1;
}

char **
backtrace_symbols(void *const *buffer, size_t size)
{
    size_t i, clen, alen, offset;
    char **rval;
    char *cp;
    Dl_info info;

    clen = size * sizeof(char *);
    rval = malloc(clen);
    if (rval == NULL)
        return NULL;

/*
    (char **)cp = &(rval[size]);
*/
/*
    cp = &(rval[size]); 
*/

    for (i = 0; i < size; i++) {
        if (dladdr(buffer[i], &info) != 0) {
            if (info.dli_sname == NULL)
                info.dli_sname = "???";
            if (info.dli_saddr == NULL)
                info.dli_saddr = buffer[i];
            offset = buffer[i] - info.dli_saddr;
            /* "0x01234567 <function+offset> at filename" */
            alen = 2 +                      /* "0x" */
                   (sizeof(void *) * 2) +   /* "01234567" */
                   2 +                      /* " <" */
                   strlen(info.dli_sname) + /* "function" */
                   1 +                      /* "+" */
                   D10(offset) +            /* "offset */
                   5 +                      /* "> at " */
                   strlen(info.dli_fname) + /* "filename" */
                   1;                       /* "\0" */
/*
            rval = realloc_safe(rval, clen + alen);
*/
            cp = malloc(alen*sizeof(char *));
            if (cp == NULL)
                return NULL;
            snprintf(cp, alen, "%p (%s+%d) at %s",
              buffer[i], info.dli_sname, offset, info.dli_fname);
        } else {
            alen = 2 +                      /* "0x" */
                   (sizeof(void *) * 2) +   /* "01234567" */
                   1;                       /* "\0" */
            rval = realloc_safe(rval, clen + alen);
            if (rval == NULL)
                return NULL;
            snprintf(cp, alen, "%p", buffer[i]);
        }
        rval[i] = cp;
        cp += alen;
    }
//  printf ("exiting !\n");
    return rval;
}

void
backtrace_symbols_fd(void *const *buffer, int size, int fd)
{
    int i, len, offset;
    char *buf;
    Dl_info info;

    for (i = 0; i < size; i++) {
        if (dladdr(buffer[i], &info) != 0) {
            if (info.dli_sname == NULL)
                info.dli_sname = "???";
            if (info.dli_saddr == NULL)
                info.dli_saddr = buffer[i];
            offset = buffer[i] - info.dli_saddr;
            /* "0x01234567 <function+offset> at filename" */
            len = 2 +                      /* "0x" */
                  (sizeof(void *) * 2) +   /* "01234567" */
                  2 +                      /* " <" */
                  strlen(info.dli_sname) + /* "function" */
                  1 +                      /* "+" */
                  D10(offset) +            /* "offset */
                  5 +                      /* "> at " */
                  strlen(info.dli_fname) + /* "filename" */
                  2;                       /* "\n\0" */
            buf = alloca(len);
            if (buf == NULL)
                return;
            snprintf(buf, len, "%p <%s+%d> at %s\n",
              buffer[i], info.dli_sname, offset, info.dli_fname);
        } else {
            len = 2 +                      /* "0x" */
                  (sizeof(void *) * 2) +   /* "01234567" */
                  2;                       /* "\n\0" */
            buf = alloca(len);
            if (buf == NULL)
                return;
            snprintf(buf, len, "%p\n", buffer[i]);
        }
        write(fd, buf, len - 1);
    }
}

#endif
