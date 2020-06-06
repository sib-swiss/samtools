
/* ==============================================================================


	Specific additions built on top of samtools to visualize several
	aligned bam files or GTL files directly in the terminal.
	by N.Guex and C.Iseli 2014-2019

    Copyright (C) SIB  - Swiss Institute of Bioinformatics,   2014-2019 Nicolas Guex and Christian Iseli
    Copyright (C) UNIL - University of Lausanne, Switzerland       2019 Nicolas Guex and Christian Iseli


    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.


	Code:       Nicolas Guex and Christian Iseli
	Contacts:   Nicolas.Guex@unil.ch and Christian.Iseli@unil.ch
	Repository: https://github.com/sib-swiss/samtools



   ============================================================================== */

/*

this source file contains portions of code taken from

  sam_view.c -- SAM<->BAM<->CRAM conversion.

	Copyright (C) 2009-2014 Genome Research Ltd.
	Portions copyright (C) 2009, 2011, 2012 Broad Institute.

	Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notices and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <inttypes.h>
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include <sys/time.h>
#include <termios.h>
#include <sys/ioctl.h>
#include <signal.h>
#include <dirent.h>
#include <pthread.h>
#include <ctype.h>
#include <sys/types.h>
#include <regex.h>

// --------------------------------------------------------------------------------------------------
//#define ROBIN  /* specific version made for Robin Liechti to support a tool in use at the HUG */
// --------------------------------------------------------------------------------------------------

#ifdef ROBIN
 #define kHalfGenomicChunk 2048  // max without changing the way I grab the genomic seq with fetch is 15000...
 #define kMaxTagsXScreens 524288
 #define CLEVER_STORE_INSERTIONS
static int gRobinPrintEveryNthLine = 1;
#else
 #define kHalfGenomicChunk 10000  // max without changing the way I grab the genomic seq with fetch is 15000...
 #define kMaxTagsXScreens 0x400000 /* guesstimate; this eats roughly 64 gigs of virtual mem (and real mem if reads are loaded) */
#endif

#define kGenomicChunk (kHalfGenomicChunk*2)
#define kMaxLineBuf 32767
#define kMaxDir 1024
#define kMaxFileName 512
#define kMaxBamName 64

#ifdef ROBIN
 #define kMaxTagLength 1024
 #define kMaxScreenCol 1024
#else
 #define kMaxTagLength 8192
 #define kMaxScreenCol 1024
#endif

#define kMaxTagInsertion 256
#define kMaxTagnameLength 48
#define kMaxTagOrdinalLength 12
#define kMaxVirutalScreen 24
#define kScreenPanelMargin 5
#define kColOffset 5
#define kMaxPatients 1024
#define kIsPairedEnd 0xFFFF
#define kIsBAM   0x1
#define kIsGTL   0x2
#define kIsTB    0x4
#define kIsMPEGG 0x8

#define kShowTagName 1
#define kShowQuality 2
#define kMaskLowQuality 4
#define kShowSNP 8 // FIXME - probably not so useful to turn this off ?  could be always on... free the v key
#define kSplitScreen 16
#define kDisplayPatientID 32
#define kZoomActivePanel 64
#define kExplorePhenotype 128
#define kHideSameMapping 256
#define kShowTagOrdinal 512
#define kFilterAllele 1024
#define kShowCoverage 2048
#define kZoomCoverage 4096
#define kReadCoverage 8192
#define kFilterClone 16384
#define kShowSequence 32768
#define kVirtScreenMaxLines 768
// --------------------------------------------------------------------------------------------------

static unsigned int gNbVirtualScreen;
static unsigned int gMaxTags;
static int gNoCompare;
static int gNoChrInBAM;
static char gGTLgenome[kMaxFileName];
static char gGTLbindir[kMaxFileName];
static char alignerConfigFile[kMaxFileName];
static char *gShowAnnotation;

typedef struct _chr_pos_t {
	unsigned int offset;
	unsigned int len;
	char name[32];
} chr_pos_t;
chr_pos_t chrPos[32];
unsigned int nbChrPos;

static char gFilterAlleleNT;
static unsigned int gFilterAllelePos;
static void usage(void);
static void help(void);
static void SIGQUITcatch (int sig);

//#define DEBUGOUTPUT
#ifdef DEBUGOUTPUT
static FILE *debug;
#endif
// --------------------------------------------------------------------------------------------------


typedef	struct PATIENT_struct  PATIENT;
struct PATIENT_struct
{
	char bamfile[32];
	char patient[16];
	int color;
};


typedef	struct TAG_struct  TAG;
struct TAG_struct
{
	unsigned long ordinal;
	char tagname[kMaxTagnameLength];
	char tagseq[kMaxTagLength];
	char tagqual[kMaxTagLength];
	char taginsertion[kMaxTagInsertion];	// FIXME: hack to store original tag
	char taginsertionQual[kMaxTagInsertion];	// FIXME:  hack to store original tag quality.
	unsigned int tagpos;
	int tagpair;   // in case of single end, a value of 1 indicate sit maps on the reverse strand
	int tagLength;
	int flag;
	int screenline;
};


typedef	struct SNPTABLE_struct  SNPTABLE;
struct SNPTABLE_struct
{
	int chr;
	int pos;
	int bgCode;
	int fgCode;
};
typedef struct _NAV_struct {
	SNPTABLE *loc;
	unsigned int max;
	unsigned int cur;
	unsigned int nb;
} NAV;
static	NAV SNPpos;
static	NAV recSNPpos;
static	NAV geneSearchPos;
static	NAV genomeSearchPos;
static	NAV *curNPpos;


#define kCoverageElements 7

typedef	struct BAMFILE_struct  BAMFILE;
struct BAMFILE_struct
{
	samFile *in;
	bam1_t *b;
	bam_hdr_t *header;
	hts_idx_t *idx;
	hts_itr_t *iter;
	unsigned short coverage[kCoverageElements * (kGenomicChunk+kMaxScreenCol)];
	unsigned int maxCoverage;
	unsigned int medianCoverage;
	unsigned int mostFrequentFragmentSize;
	unsigned int startPos;
	unsigned int retained;
	unsigned int readTags;
	int curbam;
	char patientname[16];
	char ADNIbamname[24];
	char displayname[16];
	char info[64];
	char fn[kMaxFileName];
	int cursnp;
	int color;
	int fileformat;
};

typedef	struct VIRTUALSCREEN_struct  VIRTUALSCREEN;
struct VIRTUALSCREEN_struct
{
	char *virtualscreen;
	unsigned int linelength[kVirtScreenMaxLines];
	unsigned int printed;
	unsigned int linecnt;
};

typedef	struct EXECUTIONPLAN_struct  EXECUTIONPLAN;
struct EXECUTIONPLAN_struct
{
	BAMFILE *bam;
	TAG *tag;
	int chr;
	unsigned int pos;
	int rslt;
	int screenpanel;
};




static int verbose = 0;
static char genome[kMaxLineBuf+1];
static char gAnnotationFwd[kMaxLineBuf+1];
static char gAnnotationRev[kMaxLineBuf+1];
static char gAnnotationGene[kMaxLineBuf+1];
static int volatile keeprunning = 1;
static int showmode = kShowSNP | kShowCoverage | kFilterClone;
static int gMaxTagNameLength = 0;
static int gScreenCnt = 1;
static int singleEnd = 0;
static char gLowQualityFilterThreshold = '#';
static char searchGenome[512];

#ifdef ROBIN
static int directDump = 200;
#else
static int directDump = 0;
#endif

// --------------------------------------------------------------------------------------------------
// code adapted from bam_read1 that can be found in  ../htslib/sam.c
// warning: removed the big endian test.
//
static int bam_readADNI(BGZF *fp, bam1_t *b)
{
	bam1_core_t *c = &b->core;
	int32_t block_len, ret;
	uint32_t x[8];
	if ((ret = bgzf_read(fp, &block_len, 4)) != 4)
	{
		if (ret == 0)
			return -1; // normal end-of-file
		else
			return -2; // truncated
	}
	if (bgzf_read(fp, x, 32) != 32)
	return -3;

	c->tid = x[0];
	c->pos = x[1];
	c->bin = x[2]>>16;
	c->qual = x[2]>>8&0xff;
	c->l_qname = x[2]&0xff;
	c->flag = x[3]>>16;
	c->n_cigar = x[3]&0xffff;
	c->l_qseq = x[4];
	c->mtid = x[5];
	c->mpos = x[6];
	c->isize = x[7];
	b->l_data = block_len - 32;
	if (b->l_data < 0 || c->l_qseq < 0)
		return -4;
	if ((char *)bam_get_aux(b) - (char *)b->data > b->l_data)
		return -4;
	if (b->m_data < b->l_data)
	{
		b->m_data = b->l_data;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
		if (!b->data)
			return -4;
	}
	if (bgzf_read(fp, b->data, b->l_data) != b->l_data)
		return -4;

	return 4 + block_len;

} /* bam_readADNI */
// --------------------------------------------------------------------------------------------------
static void
getCoverage(int pos, BAMFILE *bam, TAG *tag, int *lastP)
{
	unsigned int sp = 0;
	int p = tag->tagpos - pos;
	if (*lastP != 0 && p > *lastP && !(showmode & kReadCoverage))
	{
		int i;
		int max = p;
		if (max > kGenomicChunk+kMaxScreenCol)
			max = kGenomicChunk+kMaxScreenCol;
		for (i = *lastP; i < max; i++)
			bam->coverage[i * kCoverageElements + 6] += 1;
	}
	while (tag->tagseq[sp] != 0)
	{
		if (p >= kGenomicChunk+kMaxScreenCol)
			break;
		if (p >= 0)
		{
			unsigned int idx = p * kCoverageElements;
			// do not count twice when the 2 fragments overlap
			if (p >= *lastP)
			{
				unsigned int v = 0;
				switch (tag->tagseq[sp]) {
					case 'A': break;
					case 'a': break;
					case 'C': v = 1; break;
					case 'c': v = 1; break;
					case 'G': v = 2; break;
					case 'g': v = 2; break;
					case 'T': v = 3; break;
					case 't': v = 3; break;
					case '_': v = 4; break;
					case '*': v = 5; break;
					default: v = 6;
				}
				bam->coverage[idx + v] += 1;
			}
		}
		sp += 1;
		p += 1;
	}
	if (p > 0)
		*lastP = p;
	else
		*lastP = 1;
}
static int compare_short(const void *a, const void *b)
{
	short sa = *(short *)a;
	short sb = *(short *)b;
	if (sa < sb)
		return -1;
	if (sa > sb)
		return 1;
	return 0;
}
// --------------------------------------------------------------------------------------------------
// code adapted from sam_read1 that can be found in  ../htslib/sam.c
//
static int sam_readADNI(BAMFILE *bam,bam1_t *b,TAG *tag, int pos, int maxpos)
{
	int curChr = -1;		// when negative, means error.
	int r;
	unsigned int retained = 0;
	unsigned int readTags = 0;
	int dbg = 0;

	if (verbose >= 9  && dbg == 1)
		fprintf(stderr,"entering sam_readADNI maxpos = %d for %s\n",maxpos,bam->fn);

	if (bam->iter) // jump to first entry of requested chromosome.
	{
		r=sam_itr_next(bam->in, bam->iter, b);
		goto process;
	}
	while ((r=bam_readADNI((bam->in)->fp.bgzf, b)) >= 0)
	{
		process:;

		if (r >= 0)
		{
			if (b->core.tid  >= (bam->header)->n_targets || b->core.tid  < -1 ||
				b->core.mtid >= (bam->header)->n_targets || b->core.mtid < -1)
			{
				curChr = -3;
				goto done;
			}
		}

		/* ------- do our processing ---------- */

		const bam1_core_t *c = &b->core;
		if (c->qual >= 0 && c->n_cigar && c->l_qseq) // only keep good quality with cigar info and sequence
		{
			int i;

			if (c->tid != curChr) // chr change.
			{
				if (curChr != -1)
					goto done;
				curChr = c->tid;
			}
			if (verbose >= 9  && dbg == 1)
				fprintf(stderr,"readTags = %d; c->pos = %zd\n",readTags, (ssize_t) (c->pos));

#ifdef DEBUGOUTPUT
fprintf(debug,"readTags = %d; c->pos = %d\n",readTags, c->pos);
#endif

			if (c->pos > maxpos)
				goto done;

			readTags++;
			uint32_t *cigar = bam_get_cigar(b);
			for (i = 0; i < b->core.n_cigar; i++)		// do not treat I and P for now.
			{
		//		if (bam_cigar_opchr(cigar[i]) == 'I')
		//			break;
				if (bam_cigar_opchr(cigar[i]) == 'P')
					break;
			}
			if (i == b->core.n_cigar) // ok to keep.
			//if ((bam_cigar_opchr(cigar[0]) == 'M') && (bam_cigar_oplen(cigar[0]) == kTagLength))
			{
				if (b->core.l_qname > gMaxTagNameLength)	// WARN multiple threads can concurr here to alter this variable.
					gMaxTagNameLength = b->core.l_qname;
				strncpy(tag[retained].tagname,bam_get_qname(b),b->core.l_qname);
				for (i = (b->core.l_qname-1); i< kMaxTagnameLength-1/*gMaxTagNameLength*/;i++)
					tag[retained].tagname[i] = ' ';
				tag[retained].tagname[kMaxTagnameLength-1] = '\0';
				//tag[retained].tagname[b->core.l_qname+2] = '\0';
				uint8_t *s = bam_get_seq(b);
				uint8_t *q = bam_get_qual(b);
				//if (s[0] == 0xff) kputc('*', str);    case should never happen ?

			int dstp = 0;
			int srcp = 0;
			int ok = 0;
			tag[retained].ordinal = 0;
			tag[retained].taginsertion[0] = '\0';
			tag[retained].taginsertionQual[0] = '\0';
#ifdef CLEVER_STORE_INSERTIONS
			int ip=0;
#endif
			for (i = 0; i < c->n_cigar; ++i)
			{
				int k;

				switch(bam_cigar_opchr(cigar[i]))
				{

					case 'H':
						for (k = 0; k < bam_cigar_oplen(cigar[i]); k++)
						{
							srcp++;
						}
					break;
					case 'M':
						for (k = 0; k < bam_cigar_oplen(cigar[i]); k++)
						{
							tag[retained].tagseq[dstp] = "=ACMGRSVTWYHKDBN"[bam_seqi(s, srcp)];
							tag[retained].tagqual[dstp] = (char)(q[srcp] + 33);
							dstp++;
							srcp++;
							if (dstp ==	kMaxTagLength)
								goto ignore;

						}
					break;
					case 'S': // soft clipping. coordinate of match is
						for (k = 0; k < bam_cigar_oplen(cigar[i]); k++)
						{
							srcp++;
						}
					break;

					case 'D':
						for (k = 0; k < bam_cigar_oplen(cigar[i]); k++)
						{
							tag[retained].tagseq[dstp] = '_';
							tag[retained].tagqual[dstp] = '_';
							dstp++;
							if (dstp ==     kMaxTagLength)
							{
								printf("Warning, tag %s ignored expanded length kMaxTagLength\n",tag[retained].tagname);
								goto ignore;
							}
						}
					break;
					case 'N':
						for (k = 0; k < bam_cigar_oplen(cigar[i]); k++)
						{
							tag[retained].tagseq[dstp] = 'i';
							tag[retained].tagqual[dstp] = 'i';
							dstp++;
							if (dstp ==	kMaxTagLength)
							{
								printf("Warning, tag %s ignored expanded length kMaxTagLength N %d\n",tag[retained].tagname, bam_cigar_oplen(cigar[i]));
								goto ignore;
							}
						}
					break;
					case 'I':
						//printf("Warning, tag %s has insertion\n",tag[retained].tagname);
						// should find a way to store the insertion !!!
#ifdef CLEVER_STORE_INSERTIONS
						// mark presence of new insertion stored after the '|' mark.
						tag[retained].taginsertion[ip] = '|';
						tag[retained].taginsertionQual[ip] = '|';
						ip++;
#endif
						for (k = 0; k < bam_cigar_oplen(cigar[i]); k++)
						{
#ifdef CLEVER_STORE_INSERTIONS
							tag[retained].taginsertion[ip] = "=ACMGRSVTWYHKDBN"[bam_seqi(s, srcp)];
							tag[retained].taginsertionQual[ip] = (char)(q[srcp] + 33);
							ip++;
#endif
							srcp++;
						}
#ifdef CLEVER_STORE_INSERTIONS
						tag[retained].taginsertion[ip] = '\0';
						tag[retained].taginsertionQual[ip] = '\0';
#else
						// proposed way to store insertion, not particularly clever: save original tag.
						for (k=0; k < c->l_qseq;k++)
						{
							tag[retained].taginsertion[k] = "=ACMGRSVTWYHKDBN"[bam_seqi(s, k)];
							tag[retained].taginsertionQual[k] = (char)(q[k] + 33);
						}
						tag[retained].taginsertion[k] = '\0';
						tag[retained].taginsertionQual[k] = '\0';
#endif
						tag[retained].tagseq[dstp-1] = '*';
						tag[retained].tagqual[dstp-1] = '*';
					break;

					default:
						goto ignore;
					break;
				}
			}
			ok = 1;
		ignore:
				if (ok)
				{
					tag[retained].tagseq[dstp] = '\0';
					tag[retained].tagqual[dstp] = '\0';
					tag[retained].tagpos = c-> pos;
					tag[retained].tagLength = dstp; //c->l_qseq;
					tag[retained].flag = c->flag;
					if (!singleEnd)
					{
						tag[retained].tagpair = -1;
						if (c->mpos <= c->pos) // worth looking for the mate...
						{
							if (verbose >= 9  && dbg == 1)
								fprintf(stderr,"%zd <= %zd\n",(ssize_t) (c->mpos),(ssize_t) (c->pos));

							i = retained-1;
							while (i >= 0)
							{
								if (verbose >= 9  && dbg == 1)
									fprintf(stderr,"'%s'\n'%s'\n\n",tag[i].tagname,tag[retained].tagname);
								if (strcmp(tag[i].tagname,tag[retained].tagname) == 0)
								{
									tag[i].tagpair = retained;
									tag[retained].tagpair = i;
									break;
								}
								i--;
							}
						}
						if (verbose >= 9  && dbg == 1)
							fprintf(stderr,"%u\t%d\t%s\t%s\t%u\t%u\t%.20s; mpos %zd pos %zd\n",retained,tag[retained].tagpair,tag[retained].tagname,(bam->header)->target_name[c->tid],tag[retained].tagpos,c->qual,tag[retained].tagseq,(ssize_t) (c->mpos),(ssize_t) (c->pos)); 					}
					else
					{
						if (b->core.flag & BAM_FREVERSE)
							tag[retained].tagpair = 1;
						else
							tag[retained].tagpair = 0;
					}
					retained++;
				} // ok
				if (retained == gMaxTags)
				{
					fprintf(stderr,"ERROR: reached max tags (%d) ; should realloc buffer !!!!!  %d\n",retained,gMaxTags);
					fflush(stdout);
					exit(1);
					break;
				}
			}
		}
	} // while

	done:;
	memset(bam->coverage,0,kCoverageElements * (kGenomicChunk+kMaxScreenCol) * sizeof(short));
	for (r = 0; r < retained; r++)
	{
		int lastP = 0;
		if (!singleEnd)
		{
			if (tag[r].tagpair == -1)
			{
				getCoverage(pos, bam, tag + r, &lastP);
			}
			else if (tag[r].tagpair > r)
			{
				getCoverage(pos, bam, tag + r, &lastP);
				getCoverage(pos, bam, tag + tag[r].tagpair, &lastP);
			}
		}
		else
		{
			getCoverage(pos, bam, tag + r, &lastP);
		}
	}
	bam->startPos = pos;
	bam->maxCoverage = 0;
	short tem[kGenomicChunk+kMaxScreenCol];
	unsigned int tc = 0;
	for (r = 0; r < kCoverageElements * (kGenomicChunk+kMaxScreenCol); r += kCoverageElements)
	{
		unsigned int cov = bam->coverage[r] + bam->coverage[r+1] + bam->coverage[r+2] + bam->coverage[r+3] + bam->coverage[r+4] + bam->coverage[r+5] + bam->coverage[r+6];
		// put some figure.. enough to get out of noise in exome, but not too much...
		if (cov > 7)
			tem[tc++] = cov;
		if (cov > bam->maxCoverage)
			bam->maxCoverage = cov;
	}
	qsort(tem,tc,2,compare_short);
	bam->medianCoverage = 1;
	if (tc > 0)
		bam->medianCoverage = tem[tc/2];

	bam->retained = retained;
	bam->readTags = readTags;

	return (curChr);

} /* sam_readADNI */
// --------------------------------------------------------------------------------------------------
static char *bf;
static int compareFilenames(const void *a, const void *b)
{
	return strcmp(&bf[(*(int*)a) *kMaxBamName], &bf[(*(int*)b) *kMaxBamName]);
}
// --------------------------------------------------------------------------------------------------

static void sortbamNames(char *bamfiles,int bamcnt)
{
	int idx[kMaxDir];
	int i;
	char *tmp = NULL;
	int err;


	err = posix_memalign((void **)&tmp, 64, bamcnt*kMaxBamName*sizeof(char));
	if (err == 0)
	{
		for (i = 0; i< bamcnt; i++)
		{
			idx[i] = i;
			memcpy(&tmp[i*kMaxBamName],&bamfiles[i*kMaxBamName],kMaxBamName);
		}
		bf = &bamfiles[0];
		qsort(idx, bamcnt, sizeof(int), compareFilenames);

		for (i = 0; i< bamcnt; i++)
		{
			memcpy(&bamfiles[i*kMaxBamName],&tmp[idx[i]*kMaxBamName],kMaxBamName);
		}

		free(tmp);
	}
} /* sortbamNames */

// --------------------------------------------------------------------------------------------------
#define kCompStringLen (kMaxTagnameLength+16)  /* note add 16 to hold chr mapping position */
static char *tn;
static int compareTagnames(const void *a, const void *b)
{
	return strcmp(&tn[(*(int*)a) *kCompStringLen], &tn[(*(int*)b) *kCompStringLen]);
}
// --------------------------------------------------------------------------------------------------
// FIXME!   not crash free. no testing for stringlength etc...
static void compareAlignments(TAG *reftag, TAG *tag,int reftagcnt,int tagcnt)
{
	char *tmp = NULL;
	int refidx[reftagcnt];
	int idx[tagcnt];
	int tagidx;
	int err;
#ifdef DEBUGOUTPUT
	FILE *LOG = NULL;
#endif

	if (gNoCompare || !(showmode & kSplitScreen))
		return;

	if (gScreenCnt != 2)
		return;


#ifdef DEBUGOUTPUT
	LOG = fopen("/tmp/ADVIEWlog.txt","w");
#endif

	//----------- sort ref tags
	err = posix_memalign((void **)&tmp, 64, reftagcnt*kCompStringLen*sizeof(char));
	if (err != 0)
		return;

	for (tagidx = 0; tagidx < reftagcnt; tagidx++)
	{
		refidx[tagidx] = tagidx;
		//memcpy(&tmp[tagidx*kCompStringLen],&reftag[tagidx].tagname[0],kCompStringLen);
		sprintf(&tmp[tagidx*kCompStringLen],"%s%u",&reftag[tagidx].tagname[0],reftag[tagidx].tagpos);  // FIXME! potential crash possible ?
	}
	tn = &tmp[0];
	qsort(refidx, reftagcnt, sizeof(int), compareTagnames);
	free(tmp);


	// --------- sort tags
	err = posix_memalign((void **)&tmp, 64, tagcnt*kCompStringLen*sizeof(char));
	if (err != 0)
		return;

	for (tagidx = 0; tagidx < tagcnt; tagidx++)
	{
		idx[tagidx] = tagidx;
		//memcpy(&tmp[tagidx*kCompStringLen],&tag[tagidx].tagname[0],kCompStringLen);
		sprintf(&tmp[tagidx*kCompStringLen],"%s%u",&tag[tagidx].tagname[0],tag[tagidx].tagpos);  // FIXME! potential crash possible ?
	}
	tn = &tmp[0];
	qsort(idx, tagcnt, sizeof(int), compareTagnames);
	free(tmp);


	// ---------- print.
	int i = 0;
	tagidx = 0;
#ifdef DEBUGOUTPUT
	fprintf(LOG,"COMPARING TAGNAMES\n'%s'\n'%s'\n\n",&reftag[refidx[tagidx]].tagname[0],&tag[idx[i]          ].tagname[0]);
#endif
	while(1)
	{
		char *spr1 = strchr(&reftag[refidx[tagidx  ]].tagname[0], ' ');
		char *spr2 = strchr(&reftag[refidx[tagidx+1]].tagname[0], ' ');
		char *spt1 = strchr(   &tag[idx[i]          ].tagname[0], ' ');
		char *spt2 = strchr(   &tag[idx[i+1]        ].tagname[0], ' ');
		int lr1 = (spr1 != NULL) ? spr1 - reftag[refidx[tagidx]].tagname : kMaxTagnameLength;
		int lr2 = (spr2 != NULL) ? spr2 - reftag[refidx[tagidx+1]].tagname : kMaxTagnameLength;
		int lt1 = (spt1 != NULL) ? spt1 - tag[idx[i]].tagname : kMaxTagnameLength;
		int lt2 = (spt2 != NULL) ? spt2 - tag[idx[i+1]].tagname : kMaxTagnameLength;
		int lc = (lr1 >= lt1) ? lr1 : lt1;
		int lc1 = (lr1 >= lr2) ? lr1 : lr2;
		int lc2 = (lt1 >= lt2) ? lt1 : lt2;
		int cmp  = strncmp(&reftag[refidx[tagidx]].tagname[0],   &tag[idx[i]          ].tagname[0],lc);
		int cmp1 = strncmp(&reftag[refidx[tagidx]].tagname[0],&reftag[refidx[tagidx+1]].tagname[0],lc1);
		int cmp2 = strncmp(&tag[idx[i]        ].tagname[0],      &tag[idx[i+1]        ].tagname[0],lc2);

		// same tag, pair found in ref and in query.
		if (cmp == 0 && cmp1 == 0 && cmp2 == 0)
		{
#ifdef DEBUGOUTPUT
			fprintf(LOG,"%5d / %5d) %5d %s pos=%u =  (%d) %s pos=%u\n",tagidx,reftagcnt,refidx[tagidx],&reftag[refidx[tagidx]].tagname[0],reftag[refidx[tagidx]].tagpos,i,&tag[idx[i]].tagname[0],tag[idx[i]].tagpos);
			fprintf(LOG,"%5d / %5d) %5d %s pos=%u =  (%d) %s pos=%u\n\n",tagidx,reftagcnt,refidx[tagidx],&reftag[refidx[tagidx+1]].tagname[0],reftag[refidx[tagidx+1]].tagpos,i,&tag[idx[i+1]].tagname[0],tag[idx[i+1]].tagpos);
#endif
			if (reftag[refidx[tagidx]].tagpos == tag[idx[i]].tagpos)
			{
				memcpy(&reftag[refidx[tagidx]].tagname[0],"  ok  ",6);
				memcpy(&tag[idx[i]].tagname[0]           ,"  ok  ",6);
			}
			else
			{
				memcpy(&reftag[refidx[tagidx]].tagname[0]," ~ok~ ",6);
				memcpy(&tag[idx[i]].tagname[0]           ," ~ok~ ",6);
			}

			if (reftag[refidx[tagidx+1]].tagpos == tag[idx[i+1]].tagpos)
			{
				memcpy(&reftag[refidx[tagidx]].tagname[6],"  ok  ",6);
				memcpy(&tag[idx[i]].tagname[6]           ,"  ok  ",6);
			}
			else
			{
				memcpy(&reftag[refidx[tagidx]].tagname[6]," ~ok~ ",6);
				memcpy(&tag[idx[i]].tagname[6]           ," ~ok~ ",6);
			}

#ifdef DEBUGOUTPUT
			fprintf(LOG,"%5d / %5d) %5d %s pos=%u =  (%d) %s pos=%u\n",tagidx,reftagcnt,refidx[tagidx],&reftag[refidx[tagidx]].tagname[0],reftag[refidx[tagidx]].tagpos,i,&tag[idx[i]].tagname[0],tag[idx[i]].tagpos);
			fprintf(LOG,"%5d / %5d) %5d %s pos=%u =  (%d) %s pos=%u\n\n",tagidx,reftagcnt,refidx[tagidx],&reftag[refidx[tagidx+1]].tagname[0],reftag[refidx[tagidx+1]].tagpos,i,&tag[idx[i+1]].tagname[0],tag[idx[i+1]].tagpos);
#endif
			tagidx +=2;
			i +=2;
		}
		// same tag, pair found only in ref.
		else if (cmp == 0 && cmp1 == 0 && cmp2 != 0)
		{
#ifdef DEBUGOUTPUT
			fprintf(LOG,"%5d / %5d) %5d %s pos=%u =  (%d) %s pos=%u\n",tagidx,reftagcnt,refidx[tagidx],&reftag[refidx[tagidx]].tagname[0],reftag[refidx[tagidx]].tagpos,i,&tag[idx[i]].tagname[0],tag[idx[i]].tagpos);
			fprintf(LOG,"%5d / %5d) %5d %s pos=%u 1  (%d) %s pos=%u\n\n",tagidx,reftagcnt,refidx[tagidx],&reftag[refidx[tagidx+1]].tagname[0],reftag[refidx[tagidx+1]].tagpos,i,&tag[idx[i+1]].tagname[0],tag[idx[i+1]].tagpos);
#endif
			if (reftag[refidx[tagidx]].tagpos == tag[idx[i]].tagpos)
			{
				memcpy(&reftag[refidx[tagidx]].tagname[0],"  ok  ",6);
				memcpy(&tag[idx[i]].tagname[0]           ,"  ok  ",6);
			}
			else if (reftag[refidx[tagidx+1]].tagpos == tag[idx[i]].tagpos)
			{
				memcpy(&reftag[refidx[tagidx]].tagname[6],"  ok  ",6);
				memcpy(&tag[idx[i]].tagname[6]           ,"  ok  ",6);
			}
			else
			{
				memcpy(&reftag[refidx[tagidx]].tagname[0]," ~ok~ ",6);
				memcpy(&tag[idx[i]].tagname[0]           ," ~ok~ ",6);
			}

#ifdef DEBUGOUTPUT
			fprintf(LOG,"%5d / %5d) %5d %s pos=%u =  (%d) %s pos=%u\n",tagidx,reftagcnt,refidx[tagidx],&reftag[refidx[tagidx]].tagname[0],reftag[refidx[tagidx]].tagpos,i,&tag[idx[i]].tagname[0],tag[idx[i]].tagpos);
			fprintf(LOG,"%5d / %5d) %5d %s pos=%u 1  (%d) %s pos=%u\n\n",tagidx,reftagcnt,refidx[tagidx],&reftag[refidx[tagidx+1]].tagname[0],reftag[refidx[tagidx+1]].tagpos,i,&tag[idx[i+1]].tagname[0],tag[idx[i+1]].tagpos);
#endif
			tagidx +=2;
			i++;
		}
		// same tag, pair found only in query
		else if (cmp == 0 && cmp1 != 0 && cmp2 == 0)
		{
#ifdef DEBUGOUTPUT
			fprintf(LOG,"%5d / %5d) %5d %s pos=%u =  (%d) %s pos=%u\n",tagidx,reftagcnt,refidx[tagidx],&reftag[refidx[tagidx]].tagname[0],reftag[refidx[tagidx]].tagpos,i,&tag[idx[i]].tagname[0],tag[idx[i]].tagpos);
			fprintf(LOG,"%5d / %5d) %5d %s pos=%u 2  (%d) %s pos=%u\n\n",tagidx,reftagcnt,refidx[tagidx],&reftag[refidx[tagidx+1]].tagname[0],reftag[refidx[tagidx+1]].tagpos,i,&tag[idx[i+1]].tagname[0],tag[idx[i+1]].tagpos);
#endif
			if (reftag[refidx[tagidx]].tagpos == tag[idx[i]].tagpos)
			{
				memcpy(&reftag[refidx[tagidx]].tagname[0],"  ok  ",6);
				memcpy(&tag[idx[i]].tagname[0]           ,"  ok  ",6);
			}
			else if (reftag[refidx[tagidx]].tagpos == tag[idx[i+1]].tagpos)
			{
				memcpy(&reftag[refidx[tagidx]].tagname[6],"  ok  ",6);
				memcpy(&tag[idx[i]].tagname[6]           ,"  ok  ",6);
			}
			else
			{
				memcpy(&reftag[refidx[tagidx]].tagname[0]," ~ok~ ",6);
				memcpy(&tag[idx[i]].tagname[0]           ," ~ok~ ",6);
			}

#ifdef DEBUGOUTPUT
			fprintf(LOG,"%5d / %5d) %5d %s pos=%u =  (%d) %s pos=%u\n",tagidx,reftagcnt,refidx[tagidx],&reftag[refidx[tagidx]].tagname[0],reftag[refidx[tagidx]].tagpos,i,&tag[idx[i]].tagname[0],tag[idx[i]].tagpos);
			fprintf(LOG,"%5d / %5d) %5d %s pos=%u 2  (%d) %s pos=%u\n\n",tagidx,reftagcnt,refidx[tagidx],&reftag[refidx[tagidx+1]].tagname[0],reftag[refidx[tagidx+1]].tagpos,i,&tag[idx[i+1]].tagname[0],tag[idx[i+1]].tagpos);
#endif
			tagidx++;
			i += 2;
		}
		else
		{
			if (cmp < 0)
				tagidx++;
			else
				i++;
		}

/*
		// same tag ?
		if (cmp == 0)
		{
			// same tag, same pos
			if (reftag[refidx[tagidx]].tagpos == tag[idx[i]].tagpos)
			{
				// check if 2nd tag present
				int cmp2 = strcmp(&reftag[refidx[tagidx+1]].tagname[0],&tag[idx[i+1]].tagname[0]);
				if (cmp2 == 0)
				{
					if ((reftag[refidx[tagidx]].tagpos == tag[idx[i]].tagpos) && (reftag[refidx[tagidx+1]].tagpos == tag[idx[i+1]].tagpos))
					{
						memcpy(&reftag[refidx[tagidx]].tagname[0],"  ok    ok  ",6);
						memcpy(&tag[idx[i]].tagname[0]           ,"  ok    ok  ",6);
					}
					else if ((reftag[refidx[tagidx]].tagpos == tag[idx[i]].tagpos) && (reftag[refidx[tagidx+1]].tagpos != tag[idx[i+1]].tagpos))
					{
						memcpy(&reftag[refidx[tagidx]].tagname[0],"  ok   ~ok~ ",6);
						memcpy(&tag[idx[i]].tagname[0]           ,"  ok   ~ok~ ",6);
					}
					else if ((reftag[refidx[tagidx]].tagpos != tag[idx[i]].tagpos) && (reftag[refidx[tagidx+1]].tagpos == tag[idx[i+1]].tagpos))
					{
						memcpy(&reftag[refidx[tagidx]].tagname[0]," ~ok~   ok  ",6);
						memcpy(&tag[idx[i]].tagname[0]           ," ~ok~   ok  ",6);
					}
					else
					{
						memcpy(&reftag[refidx[tagidx]].tagname[0]," ~ok~  ~ok~ ",6);
						memcpy(&tag[idx[i]].tagname[0]           ," ~ok~  ~ok~ ",6);
					}
					i+=2;
					tagidx+=2;
				}
				else if (cmp2 > 0)
				{
				}
				else
				{
				}
			}
			// same tag != pos
			else
			{
				strcmp(&reftag[refidx[tagidx]].tagname[0],&tag[idx[i+1]].tagname[0]);

			}

				int overwritepos; // decide where to write the ok (left or right tag)
				int refoverwritepos; // decide where to write the ok (left or right tag)
				int overwriteidx;
				int refoverwriteidx;
				int tp = reftag[refidx[tagidx]].tagpair;

				refoverwriteidx = refidx[tagidx];
				overwriteidx = idx[i];

				if (reftag[refidx[tagidx]].tagpos < reftag[tp].tagpos)
					refoverwritepos = 0;
				else if (reftag[refidx[tagidx]].tagpos > reftag[tp].tagpos)
				{
					refoverwritepos = 6;
					refoverwriteidx = tp;
				}
				else
					refoverwritepos = 3;

				tp = tag[idx[i]].tagpair;
				if (tag[idx[i]].tagpos < tag[tp].tagpos)
				{
					overwritepos = 0;
				}
				else if (tag[idx[i]].tagpos > tag[tp].tagpos)
				{
					// check if left tag was in fact absent from ref
					if ()
					overwritepos = 6;
					overwriteidx = tp;
				}
				else
					overwritepos = 3;

				if (reftag[refidx[tagidx]].tagpos == tag[idx[i]].tagpos)
				{
					memcpy(&reftag[refoverwriteidx].tagname[refoverwritepos],"  ok  ",6);
					memcpy(&tag[overwriteidx].tagname[overwritepos],"  ok  ",6);
				}
				else
				{
					memcpy(&reftag[refoverwriteidx].tagname[refoverwritepos]," ~ok~ ",6);
					memcpy(&tag[overwriteidx].tagname[overwritepos]," ~ok~ ",6);
					fprintf(LOG,"%5d / %5d) %5d %s pos=%u ~  (%d) %s pos=%u\n",tagidx,reftagcnt,refidx[tagidx],&reftag[refidx[tagidx]].tagname[0],reftag[refidx[tagidx]].tagpos,i,&tag[idx[i]].tagname[0],tag[idx[i]].tagpos);
				}
				//printf("%5d / %5d) %5d %s == (%d) %s\n",tagidx,reftagcnt,refidx[tagidx],&reftag[refidx[tagidx]].tagname[0],i,&tag[idx[i]].tagname[0]);
				i++;
				tagidx++;
			}
			else
			{
			}
		}
		else if (cmp < 0)
		{
			fprintf(LOG,"%5d / %5d) %5d %s pos=%u <  (%d) %s pos=%u\n",tagidx,reftagcnt,refidx[tagidx],&reftag[refidx[tagidx]].tagname[0],reftag[refidx[tagidx]].tagpos,i,&tag[idx[i]].tagname[0],tag[idx[i]].tagpos);
			tagidx++;
		}
		else
		{
			fprintf(LOG,"%5d / %5d) %5d %s pos=%u >  (%d) %s pos=%u\n",tagidx,reftagcnt,refidx[tagidx],&reftag[refidx[tagidx]].tagname[0],reftag[refidx[tagidx]].tagpos,i,&tag[idx[i]].tagname[0],tag[idx[i]].tagpos);
			i++;
		}
*/
		// we want pairs, so we need at least 2 remaining
		if (tagidx + 1 >= reftagcnt)
			break;
		if (i + 1 >= tagcnt)
			break;

	}


#ifdef DEBUGOUTPUT
	fclose(LOG);
#endif


} // compareAlignments
// --------------------------------------------------------------------------------------------------
void getPatientKindLists(char *bamfiles,int bamcnt,PATIENT *patient,int *listCN,int *listMCI,int *listAD)
{
	int i;
	int cn = 0;
	int mci = 0;
	int ad = 0;
	for (i=0; i<bamcnt ; i++)
	{
		int color = 0;
		int pid = 0;

		while(patient[pid].bamfile[0] != 0)
		{
			if (strcmp(patient[pid].bamfile,&bamfiles[i*kMaxBamName]) == 0)
			{
				color = patient[pid].color;
				break;
			}
			pid++;
		}
		switch(color)
		{
			case 1: listCN[cn++] = i; break;
			case 2: listMCI[mci++] = i; break;
			case 3: listAD[ad++] = i; break;
		}
	}

	listCN[cn] = -1;
	listMCI[mci] = -1;
	listAD[ad] = -1;

} /* getPatientKindLists */
// --------------------------------------------------------------------------------------------------
void resetscreenline(TAG *tag,int retained)
{
	int i;
	for (i = 0; i<retained;i++)
	{
		tag[i].screenline = 0;
	}
}

// --------------------------------------------------------------------------------------------------

void fileappendtagfromline(TAG *tag,int retained,int chr,int line,int wait)
{
	int i = 0;
	FILE *F1 = NULL;
	FILE *F2 = NULL;

	printf("\033[H\033[J");
	fflush(stdout);

	F1 = fopen("./dbg_R1.fq","a");
	F2 = fopen("./dbg_R2.fq","a");

	if (F1 && F2)
	{
		for (i=0; i< retained; i++)
		{
			if (tag[i].screenline == line)
			{
				char revtag[1024];
				char revqual[1024];
				int k = 0;
				int kk = 0;
				int len;

				if (tag[tag[i].tagpair].taginsertion[0])  // had insertion use original tag saved for purpose.
				{
					len =strlen(tag[tag[i].tagpair].taginsertion);
					for (k=0;k<len; k++)
					{
						revqual[k] = tag[tag[i].tagpair].taginsertionQual[len-k-1];
						switch(tag[tag[i].tagpair].taginsertion[len-k-1])
						{
							case 'A': case 'a': revtag[k] = 'T';   break;
							case 'T': case 't': revtag[k] = 'A';   break;
							case 'G': case 'g': revtag[k] = 'C';   break;
							case 'C': case 'c': revtag[k] = 'G';   break;
							default: revtag[k] = tag[tag[i].tagpair].taginsertion[len-k-1]; break;
						}
						kk++;
					}
				}
				else
				{
					len =strlen(tag[tag[i].tagpair].tagseq);
					for (k=0;k<len; k++)
					{
						if (tag[tag[i].tagpair].tagseq[len-k-1] == '_')
							continue;
						revqual[kk] = tag[tag[i].tagpair].tagqual[len-k-1];
						switch(tag[tag[i].tagpair].tagseq[len-k-1])
						{
							case 'A': case 'a': revtag[kk] = 'T';   break;
							case 'T': case 't': revtag[kk] = 'A';   break;
							case 'G': case 'g': revtag[kk] = 'C';   break;
							case 'C': case 'c': revtag[kk] = 'G';   break;
							default: revtag[kk] = tag[tag[i].tagpair].tagseq[len-k-1]; break;
						}
						kk++;
					}
				}
				revtag[kk] = '\0';
				revqual[kk] = '\0';
				if (tag[i].taginsertion[0]) // had insertion use original tag saved for purpose.
				{
					int k;
					fprintf(F1,">%s(%lu)|chr%d:%d\n",tag[i].tagname,tag[i].ordinal ,chr,tag[i].tagpos);
					for (k=0;k<strlen(tag[i].taginsertion);k++) {putc(toupper(tag[i].taginsertion[k]),F1);} putc('\n',F1);
					fprintf(F1,"+\n%s\n",tag[i].taginsertionQual);
				}
				else
				{
					int k;
					fprintf(F1,">%s(%lu)|chr%d:%d\n",tag[i].tagname,tag[i].ordinal ,chr,tag[i].tagpos);
					for (k=0;k<strlen(tag[i].tagseq);k++) {if (tag[i].tagseq[k] != '_') putc(toupper(tag[i].tagseq[k]),F1);} putc('\n',F1);
					putc('+',F1); putc('\n',F1);
					for (k=0;k<strlen(tag[i].tagseq);k++) {if (tag[i].tagseq[k] != '_') putc(tag[i].tagqual[k],F1);} putc('\n',F1);
				}
				fprintf(F2,">%s(%lu)|chr%d:%d\n%s\n+\n%s\n",tag[tag[i].tagpair].tagname,tag[tag[i].tagpair].ordinal,chr,tag[tag[i].tagpair].tagpos,revtag,revqual);
				break;
			}
		}
	}

	if (F1)
		fclose(F1);
	if (F2)
		fclose(F2);

	if (i >= retained)
		printf("no tag on line %d\n",line);
	else
		printf("appended tag from line %d.\n",line);

	if (wait)
		sleep(1);

}// fileappendtagfromline
// --------------------------------------------------------------------------------------------------
void printtagfromline(TAG *tag,int retained,int chr,int line)
{
	int c;

	printf("\033[H\033[J");
	fflush(stdout);

	int i;
	for (i=0; i< retained; i++)
	{
		if (tag[i].screenline == line)
		{
			printf(">%s(%lu)|chr%d:%d\n%s\n>%s(%lu)|chr%d:%d\n%s\n",tag[i].tagname,tag[i].ordinal,chr,tag[i].tagpos,tag[i].tagseq,tag[tag[i].tagpair].tagname,tag[tag[i].tagpair].ordinal,chr,tag[tag[i].tagpair].tagpos,tag[tag[i].tagpair].tagseq);
			break;
		}
	}
	if (i >= retained)
		printf("no tag on line %d\n",line);

	read(0,&c,1);

} // printtagfromline
// --------------------------------------------------------------------------------------------------

int getpair(int genomefirstpos,int startpos, int width, TAG *tag,int retained,unsigned int *tagidx,unsigned int *fraglen, char *line)
{
	int mismatch = 0;
	int i;
	int endpos = startpos+width;

#ifdef DEBUGOUTPUT
	fprintf(debug,"genomefirstpos=%u; startpos=%u; width=%u    getting tags %d to %d\n",genomefirstpos,startpos,width,*tagidx,retained);
#endif

	for (i = *tagidx; i<retained;i++)
	{
		if ( (showmode & kHideSameMapping) && ((strncmp(&tag[i].tagname[0],"  ok    ok  ",6) == 0)))
			continue;


		int vis = 0;

		if (singleEnd || ((tag[i].tagpair != -1) && ((tag[i].tagpos < tag[ tag[i].tagpair ].tagpos) || (tag[i].tagpos == tag[ tag[i].tagpair ].tagpos && i < tag[i].tagpair))))
		{
			int leftmostpos =       tag[i].tagpos;					// start of left tag
			int tagLength = tag[i].tagLength;
			char c = '!';

/*
FIXME!   this is fine to speed up BAM files
		but it does not work for GTL, since we have several successive unsorted blocks.
		so there is no guarantee that once a tag has its leftpos out of screen (to the right)
		we will not encounter a subsequent tag that should be visible...

			if (leftmostpos > endpos) // no more hope.
			{
				*tagidx = retained;
				return(-1);
			}
*/
			// check if any part of left tag is visible
			if (((leftmostpos >= startpos) && (leftmostpos < endpos))
			|| (((leftmostpos + tagLength) > startpos) && ((leftmostpos + tagLength) < endpos))
			|| (((leftmostpos < startpos) && ((leftmostpos + tagLength) >= endpos))))
			{
				int p;
				vis = 1;
				c = '>';
				if (singleEnd && tag[i].tagpair) // hack to indicate reverse strand.
					c = '<';
				if (leftmostpos >= startpos) // start of tag visible
				{
					int k = 0;
					p = leftmostpos;
					if (showmode & kShowQuality)
					while(p < endpos)
					{
						if ( tag[i].tagseq[k] == genome[p-startpos + (startpos-genomefirstpos)] ) { line[p-startpos] = tag[i].tagqual[k]; }
						else { line[p-startpos] = tag[i].tagqual[k]; mismatch++;}
						k++;
						p++;
						if (k == tagLength) break;
					}
					else
					while(p < endpos)
					{
						if (showmode & kShowSequence)
							line[p-startpos] = tag[i].tagseq[k];
						else
						{
							if ( tag[i].tagseq[k] == genome[p-startpos + (startpos-genomefirstpos)] ) { line[p-startpos] = c; }
							else {
								if ((showmode & kMaskLowQuality) && (tag[i].tagqual[k] <= gLowQualityFilterThreshold)) {line[p-startpos] = '#';}
								else {line[p-startpos] = tag[i].tagseq[k];}
								mismatch++;
							}
						}
						k++;
						p++;
						if (k == tagLength) break;
					}
				}
				else // only end of tag visible
				{
					int k = startpos - leftmostpos;
					p = startpos;
					if (showmode & kShowQuality)
					while(p < endpos)
					{
						if ( tag[i].tagseq[k] == genome[p-startpos + (startpos-genomefirstpos)] ) { line[p-startpos] = tag[i].tagqual[k]; }
						else { line[p-startpos] = tag[i].tagqual[k]; mismatch++;}
						k++;
						p++;
						if (k == tagLength) break;
					}
					else
					while(p < endpos)
					{
						if (showmode & kShowSequence)
							line[p-startpos] = tag[i].tagseq[k];
						else
						{
							if ( tag[i].tagseq[k] == genome[p-startpos + (startpos-genomefirstpos)] ) { line[p-startpos] = c; }
							else {
								if ((showmode & kMaskLowQuality) && (tag[i].tagqual[k] <= gLowQualityFilterThreshold)) {line[p-startpos] = '#';}
								else {line[p-startpos] = tag[i].tagseq[k];}
								mismatch++;
							}
						}
						k++;
						p++;
						if (k == tagLength) break;
					}
				}
			}

			if (!singleEnd)
			{
				// check if any part of right tag is visible
				leftmostpos = tag[ tag[i].tagpair ].tagpos;	// start of right tag
				tagLength = tag[ tag[i].tagpair ].tagLength;

				if (((leftmostpos >= startpos) && (leftmostpos < endpos))
				|| (((leftmostpos + tagLength) > startpos) && ((leftmostpos + tagLength) < endpos))
				|| (((leftmostpos < startpos) && ((leftmostpos + tagLength) >= endpos))))
				{
					int p;
					vis = 1;
					if (leftmostpos >= startpos) // start of tag visible
					{
						int k = 0;
						p = leftmostpos;
						if (showmode & kShowQuality)
						while(p < endpos)
						{
							if (tag[ tag[i].tagpair ].tagseq[k] == genome[p-startpos + (startpos-genomefirstpos)] ) { line[p-startpos] = tag[ tag[i].tagpair ].tagqual[k]; }
							else { line[p-startpos] = tag[ tag[i].tagpair ].tagqual[k]; mismatch++; }

							k++;
							p++;
							if (k == tagLength) break;
						}
						else
						while(p < endpos)
						{
							if (showmode & kShowSequence)
								line[p-startpos] = tag[ tag[i].tagpair ].tagseq[k];
							else
							{
								if (tag[ tag[i].tagpair ].tagseq[k] == genome[p-startpos + (startpos-genomefirstpos)] ) { if (line[p-startpos] != '>') line[p-startpos] = '<'; else line[p-startpos] = 'X'; }
								else {
									if ((showmode & kMaskLowQuality) && (tag[ tag[i].tagpair ].tagqual[k] <= gLowQualityFilterThreshold)) {line[p-startpos] = '#';}
									else {line[p-startpos] = tag[ tag[i].tagpair ].tagseq[k];}
									mismatch++;
								}
							}

							k++;
							p++;
							if (k == tagLength) break;
						}
					}
					else // only end of tag visible
					{
						int k = startpos - leftmostpos;
						p = startpos;
						if (showmode & kShowQuality)
						while(p < endpos)
						{
							if ( tag[ tag[i].tagpair ].tagseq[k] == genome[p-startpos + (startpos-genomefirstpos)] ) { line[p-startpos] = tag[ tag[i].tagpair ].tagqual[k]; }
							else { line[p-startpos] = tag[ tag[i].tagpair ].tagqual[k]; mismatch++; }

							k++;
							p++;
							if (k == tagLength) break;
						}
						else
						while(p < endpos)
						{
							if (showmode & kShowSequence)
								line[p-startpos] = tag[ tag[i].tagpair ].tagseq[k];
							else
							{
								if ( tag[ tag[i].tagpair ].tagseq[k] == genome[p-startpos + (startpos-genomefirstpos)] ) { if (line[p-startpos] != '>') line[p-startpos] = '<'; else line[p-startpos] = 'X'; }
								else {
									if ((showmode & kMaskLowQuality) && (tag[ tag[i].tagpair ].tagqual[k] <= gLowQualityFilterThreshold)) {line[p-startpos] = '#';}
									else {line[p-startpos] = tag[ tag[i].tagpair ].tagseq[k];}
									mismatch++;
								}
							}

							k++;
							p++;
							if (k == tagLength) break;
						}
					}
				} // right tag
			}
#ifdef DEBUGOUTPUT
			fprintf(debug,"%.50s %8u %8u %c  %c visible=%d\n",tag[i].tagname,tag[i].tagpos,tag[ tag[i].tagpair ].tagpos,(tag[i].tagpos > tag[ tag[i].tagpair ].tagpos) ? '*' : ' ' ,c,vis);
#endif
			if (showmode & kFilterAllele)
			{
				int lp = gFilterAllelePos - startpos;
				fflush(stdout);
				if ((lp > 0) && (lp < kMaxScreenCol))
				{
					if (line[lp] != gFilterAlleleNT)
					{
						memset(line,'.',kMaxScreenCol);
						vis = 0;
					}
				}
			}
		}
		else
		{
#ifdef DEBUGOUTPUT
			fprintf(debug,"SKIPPED %.50s %8u %8u\n",tag[i].tagname,tag[i].tagpos,tag[ tag[i].tagpair ].tagpos);
#endif
		}
		if (vis)
		{
			unsigned int fragend = tag[ tag[i].tagpair ].tagpos + tag[ tag[i].tagpair ].tagLength;
			*fraglen = ((fragend - tag[i].tagpos));
			break;
		}
	} // i
#ifdef DEBUGOUTPUT
	fflush(debug);
#endif

	*tagidx = i;
	return(mismatch);

} /* getpair */
// --------------------------------------------------------------------------------------------------
void dumpVirtualScreenToFile(VIRTUALSCREEN *screenpanel,FILE *f,char *name)
{
	int i;
	int from = 0;

	if (f)
	{
		for (i = 0; i < screenpanel->linecnt; i++)
		{
			fprintf(f,"%.2s%.2s\t%s\n",&screenpanel->virtualscreen[from],&screenpanel->virtualscreen[from+screenpanel->linelength[i]-2],name);
			from += screenpanel->linelength[i];
			fflush(f);
		}
	}

} /* dumpVirtualScreenToFile */
// --------------------------------------------------------------------------------------------------
void printchunkBW(TAG *tag, int retained, int chr, int line,unsigned int firstGenomeNTavail,unsigned int startpos,int width,int height,VIRTUALSCREEN *screenpanel)
{

	unsigned int tagidx;
	unsigned int printed = 0;
	char white[kMaxScreenCol];
	int ll = 1;
	int bufpos = 0;
	char *buffer = screenpanel->virtualscreen;
	unsigned int *linelength = screenpanel->linelength;
	unsigned int maxTagNameLength = gMaxTagNameLength;
	unsigned int tagNameStartPos = 0;

	int mismatchallowed = 999999;

	if (maxTagNameLength > width)
	{
		tagNameStartPos = maxTagNameLength - width;
		maxTagNameLength = width;
	}

	// print tags
	tagidx = 0;  // note this is inefficient. the whole tag array is traversed from the beginning until the first tag to print is detected...
				// we could remember print from last time (per screen) or bracket to identify start...

	while(printed < height)
	{
		unsigned int fraglen;
		memset(white,'.',kMaxScreenCol);
		if (getpair(firstGenomeNTavail,startpos, width,tag, retained,&tagidx,&fraglen,white) <= mismatchallowed)
		{
			if (tagidx < retained)
			{
				if (ll++ >= line) // print only once we skipped ll lines due to scrolling instructions.
				{
					if (!(showmode & kShowQuality))
					{
						linelength[printed] = sprintf(&buffer[bufpos],"%.*s",width,white); bufpos += linelength[printed++];
					}
					else
					{
						if (showmode & kShowTagName)
						{
							snprintf(white,maxTagNameLength,"%s",&tag[tagidx].tagname[tagNameStartPos]);
							white[maxTagNameLength-1] = ' ';
						}
						else if (showmode & kShowTagOrdinal)
						{
							snprintf(white,kMaxTagOrdinalLength,"%-*lu",kMaxTagOrdinalLength-1,tag[tagidx].ordinal);
							white[kMaxTagOrdinalLength-1] = ' ';
						}
						linelength[printed] = sprintf(&buffer[bufpos],"%4u %.*s",ll-1,width,white); bufpos += linelength[printed++];
					}
				}
			}
		}

		tagidx++;
		if (tagidx >= retained)
			break;
	}

	screenpanel->linecnt = printed;

} /* printchunkBW */
// --------------------------------------------------------------------------------------------------

void printchunkForRobin(TAG *tag, int retained, int chr,unsigned int firstGenomeNTavail,unsigned int startpos,int width)
{

	unsigned int tagidx;
	char white[kMaxScreenCol];
	int ll = 1;

	unsigned int maxTagNameLength = gMaxTagNameLength;
	unsigned int tagNameStartPos = 0;

	int mismatchallowed = 99999999;
	printf("%d",width);
	if (maxTagNameLength > width)
	{
		tagNameStartPos = maxTagNameLength - width;
		maxTagNameLength = width;
	}

	// print tags
	tagidx = 0;  // note this is inefficient. the whole tag array is traversed from the beginning until the first tag to print is detected...
				// we could remember print from last time (per screen) or bracket to identify start...

	printf("CHR%d:%u\n",chr,startpos+1); // +1 to be one base for the genome....
	printf("%.*s\n",width,&genome[startpos-firstGenomeNTavail]);
#ifdef ROBIN
int RobinOutputedLine = 0;
#endif
	while(1)
	{
		unsigned int fraglen;
		memset(white,'.',kMaxScreenCol);
		if (getpair(firstGenomeNTavail,startpos, width,tag, retained,&tagidx,&fraglen,white) <= mismatchallowed)
		{
#ifdef ROBIN
			if (++RobinOutputedLine < gRobinPrintEveryNthLine)
				goto RobinSkip;
			else
				RobinOutputedLine = 0;
#endif
			if (tagidx < retained)
			{
				if (showmode & kShowTagName)
				{
					snprintf(white,maxTagNameLength,"%s",&tag[tagidx].tagname[tagNameStartPos]);
					white[maxTagNameLength-1] = ' ';
				}
				else if (showmode & kShowTagOrdinal)
				{
					snprintf(white,kMaxTagOrdinalLength,"%-*lu",kMaxTagOrdinalLength-1,tag[tagidx].ordinal);
					white[kMaxTagOrdinalLength-1] = ' ';
				}
#ifdef CLEVER_STORE_INSERTIONS
				if (!(showmode & kShowQuality))
					printf("%.*s%s\n",width,white,&tag[tagidx].taginsertion[0]);
				else
					printf("%4u %.*s%s\n",ll-1,width,white,&tag[tagidx].taginsertionQual[0]);
#else
				if (!(showmode & kShowQuality))
					printf("%.*s\n",width,white);
				else
					printf("%4u %.*s\n",ll-1,width,white);
#endif
			}
		}
#ifdef ROBIN
RobinSkip:;
#endif

		tagidx++;
#ifdef ROBIN
		if (tagidx >= retained*gRobinPrintEveryNthLine)
			break;
#endif
	}

} /* printchunkForRobin */
// --------------------------------------------------------------------------------------------------



void flushscreen(VIRTUALSCREEN *screenpanel,int panelcnt)
{
	int i;
	unsigned int from[gNbVirtualScreen];
	for (i = 0; i< panelcnt; i++)
		from[i] = 0;

	puts("\033[f");// move to top.

	for (i = 0; i < screenpanel->linecnt; i++)
	{
		int k;
		printf("\033[%d;1f",i+1);
		for (k=0;k<panelcnt;k++)
		{
			printf("%.*s",screenpanel[k].linelength[i],&screenpanel[k].virtualscreen[from[k]]);
			from[k] += screenpanel[k].linelength[i];
		}
	}

} /* flushscreen */
// --------------------------------------------------------------------------------------------------
void flushscreenPanelToFile(char *fn, VIRTUALSCREEN *screenpanel,int panel,int lcnt)
{
	int i;
	unsigned int from = 0;

	FILE *of = fopen(fn,"w");

	for (i = 0; i < lcnt; i++)
	{
		fprintf(of,"%.*s\n",screenpanel[panel].linelength[i],&screenpanel[panel].virtualscreen[from]);
		from += screenpanel[panel].linelength[i];
	}
	fclose(of);
	
} /* flushscreenPanelToFile */
// --------------------------------------------------------------------------------------------------

void printchunk(TAG *tag, BAMFILE *bam, int chr, int line,unsigned int firstGenomeNTavail,unsigned int startpos,int width,int height,VIRTUALSCREEN *screenpanel,int color)
{
	int retained = bam->retained;
	unsigned int tagidx;
	unsigned int printed = 0;
	char white[kMaxScreenCol];
	char initBGcolor[16];
	char endBGcolor[8];
	int k;
	int i;
	int ll = 1;
	int roundstartpos = startpos;
	int bufpos = 0;
	char *buffer = screenpanel->virtualscreen;
	unsigned int *linelength = screenpanel->linelength;
	unsigned int maxTagNameLength = gMaxTagNameLength;
	unsigned int tagNameStartPos = 0;


	switch(color)
	{
		case 0: initBGcolor[0] = 0;  endBGcolor[0] = 0;  break;
		case 1: strcpy(initBGcolor,"\033[1;42m\033[1;37m");  strcpy(endBGcolor,"\033[0m");  break;  // green
		case 2: strcpy(initBGcolor,"\033[1;43m\033[1;37m");  strcpy(endBGcolor,"\033[0m");  break;  // yellow
		case 3: strcpy(initBGcolor,"\033[1;41m\033[1;37m");  strcpy(endBGcolor,"\033[0m");  break;  // red
	}

	int mismatchallowed = 999999;

	if (maxTagNameLength > width)
	{
		tagNameStartPos = maxTagNameLength - width;
		maxTagNameLength = width;
	}
	if (startpos - bam->startPos <= kGenomicChunk && (showmode & kShowCoverage) && !(showmode & kZoomCoverage))
	{
		const unsigned int cpart = (height > 34) ? 9 : 3;
		const unsigned int step = ((bam->medianCoverage * 2) > (cpart * 8)) ? (bam->medianCoverage * 2) / (cpart * 8) : 1;
		const unsigned int halfStep = (step > 1) ? step >> 1 : 1;
		const char *ntcode = "ACGT_*N";
		char *panel = malloc(cpart * width * 16 * sizeof(char));
		unsigned int *panelll = calloc(cpart,sizeof(unsigned int));
		for (i = 0; i < width; i++)
		{
			unsigned int p = startpos - bam->startPos + i;
			unsigned int cp = kCoverageElements * p;
			unsigned int cov = 0;
			unsigned int pcov = 0;
			struct { char c; unsigned int cov; } cs[7];
			unsigned int csn = 0;
			unsigned int j;
			for (j = 0; j < 6; j++)
			{
				unsigned int k = bam->coverage[cp + j];
				cov += k;
				if (k >= halfStep)
				{
					unsigned int csi = csn;
					while (csi > 0 && k > cs[csi - 1].cov)
					{
						cs[csi] = cs[csi - 1];
						csi -= 1;
					}
					csn += 1;
					cs[csi].c = ntcode[j];
					if (cs[csi].c == genome[startpos-firstGenomeNTavail+i])
						cs[csi].c = ' ';
					cs[csi].cov = k;
					pcov += k;
				}
			}
			int printed = 0;
			if (pcov > 0)
			{
				int max = ((bam->medianCoverage * 2) > (cpart * 8)) ? (bam->medianCoverage * 2) : cpart * 8;
				if (cov > max)
					max = cov;
				int lv = max / (cpart * 8);
				if (lv < 1)
					lv = 1;
				if (lv * cpart * 8 < max)
					lv += 1;
				max = lv * cpart * 8;
				j = 0;
				int sum = cs[0].cov;
				//for (j = 0; j < csn; j++)
				while (j < csn)
				{
					while (printed * lv * 8 < sum)
					{
						if (!(printed < cpart))
						{
							fprintf(stderr,"printed grows too large cpart:%d printed:%d cov:%d max:%d lv:%d sum:%d (printed * lv * 8):%d\n",cpart,printed,cov,max,lv,sum,printed * lv * 8);
							break;
						}
						int diff = (sum - printed * lv * 8) / lv;
						if (diff <= 0)
							break;
						if (diff > 8)
							diff = 8;
						char c1 = cs[j].c;
						char c2 = cs[j].c;
						if (diff < 8)
						{
							j += 1;
							if (j < csn)
							{
								sum += cs[j].cov;
								int diff2 = (sum - printed * lv * 8) / lv;
								if (diff2 >= 8)
								{
									// we add a second color...
									c2 = cs[j].c;
								}
							}
						}
						if (c1 != ' ' || c1 != c2)
						{
							char c = '0';
							unsigned int cnt = 0;
							switch(c1)
							{
								case ' ': c = '7';   break; // white
								case 'A': c = '6';   break; // cyan
								case 'T': c = '3';   break; // orange
								case 'G': c = '5';   break;
								case 'C': c = '1';   break; // red
								case '_': c = '4';   break; // blue
								case '*': c = '2';   break; // green
							}
							panel[printed * width * 16 + panelll[printed] + cnt++] = '\033';
							panel[printed * width * 16 + panelll[printed] + cnt++] = '[';
							panel[printed * width * 16 + panelll[printed] + cnt++] = '1';
							panel[printed * width * 16 + panelll[printed] + cnt++] = ';';
							panel[printed * width * 16 + panelll[printed] + cnt++] = '3';
							panel[printed * width * 16 + panelll[printed] + cnt++] = c;
							if (c1 != c2)
							{
								switch(c2)
								{
									case ' ': c = '7';   break; // white
									case 'A': c = '6';   break; // cyan
									case 'T': c = '3';   break; // orange
									case 'G': c = '5';   break;
									case 'C': c = '1';   break; // red
									case '_': c = '4';   break; // blue
									case '*': c = '2';   break; // green
								}
								panel[printed * width * 16 + panelll[printed] + cnt++] = ';';
								panel[printed * width * 16 + panelll[printed] + cnt++] = '1';
								panel[printed * width * 16 + panelll[printed] + cnt++] = '0';
								panel[printed * width * 16 + panelll[printed] + cnt++] = c;
							}
							panel[printed * width * 16 + panelll[printed] + cnt++] = 'm';
							panelll[printed] += cnt;
						}
						panel[printed * width * 16 + panelll[printed]] = '\xe2';
						panel[printed * width * 16 + panelll[printed] + 1] = '\x96';
						panel[printed * width * 16 + panelll[printed] + 2] = '\x80' + diff;
						panelll[printed] += 3;
						if (c1 != ' ' || c1 != c2)
						{
							panel[printed * width * 16 + panelll[printed]] = '\033';
							panel[printed * width * 16 + panelll[printed] + 1] = '[';
							panel[printed * width * 16 + panelll[printed] + 2] = '0';
							panel[printed * width * 16 + panelll[printed] + 3] = 'm';
							panelll[printed] += 4;
						}
						printed += 1;
					}
					j += 1;
					if (j < csn)
						sum += cs[j].cov;
				}
			}
			while (printed < cpart)
			{
				panel[printed * width * 16 + panelll[printed]] = ' ';
				panelll[printed] += 1;
				printed += 1;
			}
		}
		unsigned int val = ((bam->medianCoverage * 2) > (cpart * 8)) ? (bam->medianCoverage * 2) : cpart * 8;
		if (val > 999)
			val = 999;
		i = cpart - 1;
		linelength[printed] = sprintf(&buffer[bufpos]," %3u %.*s",val,panelll[i],panel + i * width * 16); bufpos += linelength[printed++];
		i -= 1;
		for ( ; i > (cpart / 2); i--)
		{
			linelength[printed] = sprintf(&buffer[bufpos],"     %.*s",panelll[i],panel + i * width * 16); bufpos += linelength[printed++];
		}
		val = ((bam->medianCoverage * 2) > (cpart * 8)) ? bam->medianCoverage : cpart * 4;
		if (val > 999)
			val = 999;
		linelength[printed] = sprintf(&buffer[bufpos]," %3u %.*s",val,panelll[i],panel + i * width * 16); bufpos += linelength[printed++];
		i -= 1;
		for ( ; i > 0; i--)
		{
			linelength[printed] = sprintf(&buffer[bufpos],"     %.*s",panelll[i],panel + i * width * 16); bufpos += linelength[printed++];
		}
		val = step;
		if (val > 999)
			val = 999;
		linelength[printed] = sprintf(&buffer[bufpos]," %3u %.*s",val,panelll[i],panel + i * width * 16); bufpos += linelength[printed++];
		free(panel);
		free(panelll);
	}
	// print scale
	k = 0;
	i= startpos%10;
	while(i++ < 10) { white[k++] = ' '; roundstartpos++; }
	for (i = 0; i< (width / 9); i+=2)
		k+=sprintf(&white[k],"%10u          ",roundstartpos+(i+1)*10);
	linelength[printed] = sprintf(&buffer[bufpos],"%s C%2d %.*s%s",initBGcolor,chr,width,white,endBGcolor); bufpos += linelength[printed++];

	// print tickmarks
	k = 0;
	i=startpos%10;
	while(i++ < 10) { white[k++] = ' '; }
	for (i = 0; i< (width / 9); i+=2)
		k+=sprintf(&white[k],"    .    |    .    |");
	linelength[printed] = sprintf(&buffer[bufpos],"     %.*s",width,white); bufpos += linelength[printed++];

	// print sequence
	if (searchGenome[0] != 0)
	{
		size_t len = strlen(searchGenome);
		char gc[width + len * 2 - 1];
		strncpy(gc,genome + startpos - firstGenomeNTavail - len + 1, width + len * 2 - 2);
		gc[width + len * 2 - 2] = 0;
		char *gs = gc + len - 1;
		char *s = gc;
		regex_t re;
		regmatch_t rm;
		regcomp(&re, searchGenome, REG_EXTENDED);
		if (regexec(&re, s, 1, &rm, REG_NOTBOL|REG_NOTEOL) == 0)
			s += rm.rm_so;
		else
			s = NULL;
		int cur = 0;
		char cs[width*12+1];
		cs[0] = 0;
		while (s != NULL && s < gs + width)
		{
			int covlen = rm.rm_eo - rm.rm_so;
			if (s < gs + cur)
			{
				covlen -= gs + cur - s;
				s = gs + cur;
			}
			if (s > gs + cur)
			{
				int delta = s - gs - cur;
				strncat(cs,gs + cur,delta);
				cur += delta;
			}
			strcat(cs,"\033[1;45m");
			if (covlen < width - cur)
			{
				strncat(cs,gs + cur,covlen);
				cur += covlen;
				s += covlen;
				if (regexec(&re, s, 1, &rm, REG_NOTBOL|REG_NOTEOL) == 0)
					s += rm.rm_so;
				else
					s = NULL;
			}
			else
			{
				strncat(cs,gs + cur,width - cur);
				cur = width;
				s = NULL;
			}
			strcat(cs,"\033[0m");
		}
		regfree(&re);
		if (cur < width)
			strncat(cs,gs + cur,width - cur);
		linelength[printed] = sprintf(&buffer[bufpos],"     %s",cs); bufpos += linelength[printed++];
	}
	else
	{
		linelength[printed] = sprintf(&buffer[bufpos],"     %.*s",width,&genome[startpos-firstGenomeNTavail]); bufpos += linelength[printed++];
	}
	if (gShowAnnotation != NULL)
	{
		linelength[printed] = sprintf(&buffer[bufpos],"     %.*s",width,&gAnnotationRev[startpos-firstGenomeNTavail]); bufpos += linelength[printed++];
		linelength[printed] = sprintf(&buffer[bufpos],"     %.*s",width,&gAnnotationFwd[startpos-firstGenomeNTavail]); bufpos += linelength[printed++];
		linelength[printed] = sprintf(&buffer[bufpos],"     %.*s",width,&gAnnotationGene[startpos-firstGenomeNTavail]); bufpos += linelength[printed++];
	}
	if (startpos - bam->startPos <= kGenomicChunk && (showmode & kShowCoverage) && (showmode & kZoomCoverage))
	{
		const unsigned int cpart = height - printed;
		const unsigned int step = ((bam->medianCoverage * 2) > cpart) ? (bam->medianCoverage * 2) / cpart : 1;
		const unsigned int halfStep = (step > 1) ? step >> 1 : 1;
		const char *ntcode = (showmode & kReadCoverage) ? "ACGT_*N" : "ACGT_*'";
		char panel[cpart * width];
		memset(panel,' ',cpart * width);
		for (i = 0; i < width; i++)
		{
			unsigned int p = startpos - bam->startPos + i;
			unsigned int cp = kCoverageElements * p;
			unsigned int cov = 0;
			unsigned int pcov = 0;
			struct { char c; unsigned int cov; } cs[7];
			unsigned int csn = 0;
			unsigned int j;
			for (j = 0; j < 7; j++)
			{
				unsigned int k = bam->coverage[cp + j];
				cov += k;
				if (k >= halfStep)
				{
					unsigned int csi = csn;
					while (csi > 0 && k > cs[csi - 1].cov)
					{
						cs[csi] = cs[csi - 1];
						csi -= 1;
					}
					csn += 1;
					cs[csi].c = ntcode[j];
					if (cs[csi].c == genome[startpos-firstGenomeNTavail+i])
						cs[csi].c = '`';
					cs[csi].cov = k;
					pcov += k;
				}
			}
			if (pcov == 0)
				continue;
			unsigned int max = (cov > bam->medianCoverage * 2) ? cov : bam->medianCoverage * 2;
			unsigned int lv = max / cpart;
			if (lv < 1)
				lv = 1;
			if (lv * cpart < max)
				lv += 1;
			max = lv * cpart;
			unsigned int sum = 0;
			unsigned int printed = 0;
			for (j = 0; j < csn; j++)
			{
				sum += cs[j].cov;
				while (printed * lv < sum)
				{
					if (!(printed < cpart))
					{
						fprintf(stderr,"printed grows too large cpart:%d printed:%d lv:%d sum:%d (printed * lv):%d\n",cpart,printed,lv,sum,printed * lv);
						break;
					}
					panel[printed * width + i] = cs[j].c;
					printed += 1;
				}
			}
		}
		unsigned int val = step;
		if (val > 999)
			val = 999;
		i = 0;
		linelength[printed] = sprintf(&buffer[bufpos]," %3u %.*s",val,width,panel + i * width); bufpos += linelength[printed++];
		i += 1;
		for ( ; i < (cpart / 2); i++)
		{
			linelength[printed] = sprintf(&buffer[bufpos],"     %.*s",width,panel + i * width); bufpos += linelength[printed++];
		}
		val = step * cpart / 2;
		if (val > 999)
			val = 999;
		linelength[printed] = sprintf(&buffer[bufpos]," %3u %.*s",val,width,panel + i * width); bufpos += linelength[printed++];
		i += 1;
		for ( ; i < (cpart - 1); i++)
		{
			linelength[printed] = sprintf(&buffer[bufpos],"     %.*s",width,panel + i * width); bufpos += linelength[printed++];
		}
		val = step * cpart;
		if (val > 999)
			val = 999;
		linelength[printed] = sprintf(&buffer[bufpos]," %3u %.*s",val,width,panel + i * width); bufpos += linelength[printed++];
	}

	// print tags
	tagidx = 0;  // note this is inefficient. the whole tag array is traversed from the beginning until the first tag to print is detected...
				// we could remember print from last time (per screen) or bracket to identify start...

	resetscreenline(tag,retained);
	while(printed < height)
	{
		unsigned int fraglen;
		memset(white,'.',kMaxScreenCol);
		if (getpair(firstGenomeNTavail,startpos, width,tag, retained,&tagidx,&fraglen,white) <= mismatchallowed)
		{
			if (tagidx < retained)
			{
				if (ll++ >= line) // print only once we skipped ll lines due to scrolling instructions.
				{
					char s[12*kMaxScreenCol];	// some extra space for escape codes...

					if (showmode & kShowTagName)
					{
						snprintf(s,maxTagNameLength,"%s",&tag[tagidx].tagname[tagNameStartPos]);
						s[maxTagNameLength-1] = ' ';
						k = maxTagNameLength;
						i = maxTagNameLength;
					}
					else if (showmode & kShowTagOrdinal)
					{
						snprintf(s,kMaxTagOrdinalLength,"%-*lu",kMaxTagOrdinalLength-1,tag[tagidx].ordinal);
						s[kMaxTagOrdinalLength-1] = ' ';
						k = kMaxTagOrdinalLength;
						i = kMaxTagOrdinalLength;
					}
					else
					{
						i = 0;
						k =0;
					}

					// colorize in function of fraglen
					char lineNum[32];
					if (fraglen >= bam->mostFrequentFragmentSize - bam->mostFrequentFragmentSize/2)
					{
						if (fraglen >= bam->mostFrequentFragmentSize + bam->mostFrequentFragmentSize/2)
							sprintf(lineNum,"\033[1;36m%4d\033[0m",ll-1);
						else
							sprintf(lineNum,"%4d",ll-1);
					}
					else
						sprintf(lineNum,"\033[1;31m%4d\033[0m",ll-1);  // red


					if (!(showmode & kShowQuality))
					{
						for (; i<width;i++)
						{
							char c = white[i];
							switch(c)
							{
								case 'A': k += sprintf(&s[k],"\033[1;36m%c\033[0m",c);   break; // cyan
								case 'T': k += sprintf(&s[k],"\033[1;33m%c\033[0m",c);   break; // orange
								case 'G': k += sprintf(&s[k],"\033[1;35m%c\033[0m",c);   break;
								case 'C': k += sprintf(&s[k],"\033[1;31m%c\033[0m",c);   break; // red
								case '_': k += sprintf(&s[k],"\033[1;47m%c\033[0m",c);   break; // white
								case '*': k += sprintf(&s[k],"\033[1;47m%c\033[0m",c);   break; // white
						//		case '>': s[k++] = '>';   break;
						//		case '<': s[k++] = '<';   break;
								default: s[k++] = c; break;
							}
						}
						linelength[printed] = sprintf(&buffer[bufpos],"%s %.*s",lineNum,k,s); bufpos += linelength[printed++];
					}
					else
					{
						if (showmode & kShowTagName)
						{
							snprintf(white,maxTagNameLength,"%s",&tag[tagidx].tagname[tagNameStartPos]);
							white[maxTagNameLength-1] = ' ';
						}
						else if (showmode & kShowTagOrdinal)
						{
							snprintf(white,kMaxTagOrdinalLength,"%-*lu",kMaxTagOrdinalLength-1,tag[tagidx].ordinal);
							white[kMaxTagOrdinalLength-1] = ' ';
						}
						linelength[printed] = sprintf(&buffer[bufpos],"%s %.*s",lineNum,width,white); bufpos += linelength[printed++];
					}
					tag[tagidx].screenline = ll-1;
				}
			}
		}

		tagidx++;
		if (tagidx >= retained)
			break;
	}
	memset(white,'.',kMaxScreenCol);
	screenpanel->printed = printed;
	while(printed < height)
	{
		ll++;
		linelength[printed] = sprintf(&buffer[bufpos],"%4u %.*s",ll-1,width,white); bufpos += linelength[printed++];
	}
	screenpanel->linecnt = printed;

} /* printchunk */
// --------------------------------------------------------------------------------------------------


static int openBAM(BAMFILE *bam, PATIENT *p)
{
	int pid;


	/* ---- check type from name --- */
#ifndef ROBIN
	int len = strlen(bam->fn);
	if ((len < 4) || (strncmp(&bam->fn[len-4],".bam",4) != 0))
	{

		if (strncmp(&bam->fn[len-3],".rl",3) == 0)
		{
			bam->fileformat = kIsTB;
		}
		else if (strncmp(&bam->fn[len-4],".mgg",4) == 0)
		{
			bam->fileformat = kIsMPEGG;
		}
		else
		{
			// assume GTL
			bam->fileformat = kIsGTL;
		}

		int from = len-16;
		if (from < 0)
			from = 0;

		memcpy(bam->displayname,&bam->fn[from],16);

		pid = 0;
		bam->patientname[0] = 0;
		while(p[pid].bamfile[0] != 0)
		{
			if (strcmp(bam->ADNIbamname,p[pid].bamfile) == 0)
			{
				memcpy(bam->patientname,p[pid].patient,12);
				bam->color = p[pid].color;
				break;
			}
			pid++;
		}
		return(0);
	}
#endif
 	/* --------- open BAM file ----------------*/

	if ((bam->in = sam_open(bam->fn, "r")) == 0)
	{
		fprintf(stderr, "failed to open \"%s\" for reading\n", bam->fn);
		return(1);
	}

	if ((bam->header = sam_hdr_read(bam->in)) == 0)
	{
		fprintf(stderr, "[ADVIEW] failed to read the header from \"%s\".\n", bam->fn);
		return(1);
	}

	if ((bam->in)->is_be)
	{
		printf("ERROR, big endian not supported\n");
		return(1);
	}

	/* --------- open index ----------------*/

	bam->idx = sam_index_load(bam->in, bam->fn); // load index
	if (bam->idx == 0)
	{
		fprintf(stderr, "index of %s not found\n",bam->fn);
		return(1);
	}

	int l = strlen(bam->fn);
	int from = l-16;
	if (from < 0)
		from = 0;

	strncpy(bam->displayname,&bam->fn[from],12);

	if (p[0].bamfile[0])
	{
		pid = 0;
		bam->patientname[0] = 0;
		while(p[pid].bamfile[0] != 0)
		{
			if (strcmp(bam->ADNIbamname,p[pid].bamfile) == 0)
			{
				int from = strlen(p[pid].patient)-16;
				if (from < 0)
					from = 0;

				memcpy(bam->patientname,&p[pid].patient[from],12);
				bam->color = p[pid].color;
				break;
			}
			pid++;
		}
	}
	bam->fileformat = kIsBAM;
	return(0);

} /* openBAM */
// --------------------------------------------------------------------------------------------------
static void closeBAM(BAMFILE *bam)
{
	if (bam->idx)
		hts_idx_destroy(bam->idx); // destroy the BAM index

	if (bam->in)
		sam_close(bam->in);

	if (bam->header)
		bam_hdr_destroy(bam->header);

	bam->in =NULL;
	bam->idx =NULL;
	bam->header =NULL;
	bam->color = 0;
	bam->fileformat = 0;

} /* closeBAM */
// --------------------------------------------------------------------------------------------------

static int getGTLchunk(char *fn,int chr,int pos,BAMFILE *bam,TAG *tag,int realign)
{


	FILE *f[4];
	char cmd[kMaxFileName*4];
	char buf[8192];
	int i,t;
	char tagKind[4] = {'p','n','m','a'};
	unsigned int fragSizeCnt[1024];
	unsigned int prev1 = 0, prev2 = 0;
	unsigned int fCnt;
	
	
	i = 0;
	memset(&fragSizeCnt[0],0,1024);
	memset(bam->coverage,0,kCoverageElements * (kGenomicChunk+kMaxScreenCol) * sizeof(short));

	// open files
	if (bam->fileformat == kIsGTL)
	{
		fCnt = 4;
		for (t = 0; t<fCnt;t++)
		{
			if (realign && ((tagKind[t] == 'a') || (tagKind[t] == 'm')))
			{
				if (gGTLbindir[0] == 0)
					sprintf(cmd,"GTLdecompress -g %s -i %s -C %d -P %d..%d -%c -o realign -R %s",gGTLgenome,fn,chr,pos-1250,pos+1250+kMaxScreenCol,tagKind[t],alignerConfigFile);
				else
					sprintf(cmd,"%s/GTLdecompress -g %s -i %s -C %d -P %d..%d -%c -o realign -R %s",gGTLbindir,gGTLgenome,fn,chr,pos-1250,pos+1250+kMaxScreenCol,tagKind[t],alignerConfigFile);
			}
			else
			{
				if (gGTLbindir[0] == 0)
					sprintf(cmd,"GTLdecompress -g %s -i %s -C %d -P %d..%d -%c -o ADNIview",gGTLgenome,fn,chr,pos,pos+kGenomicChunk+kMaxScreenCol,tagKind[t]);
				else
					sprintf(cmd,"%s/GTLdecompress -g %s -i %s -C %d -P %d..%d -%c -o ADNIview",gGTLbindir,gGTLgenome,fn,chr,pos,pos+kGenomicChunk+kMaxScreenCol,tagKind[t]);
			}
			//sleep(2);
			if (verbose)
				fprintf(stderr, "%s\n",cmd);
			f[t] = popen (cmd, "r");
			if (!f[t])
			{
				fprintf(stderr, "Error:Cannot execute command %s\n",cmd);
				return(-1);
			}
		}
	}
	else if (bam->fileformat == kIsMPEGG)
	{
		fCnt = 4;
		for (t = 0; t<fCnt;t++)
		{
			if (gGTLbindir[0] == 0)
				sprintf(cmd,"MPEGGdecompress -i %s -A chr%d:%d-%d -%c",fn,chr,pos,pos+kGenomicChunk+kMaxScreenCol,tagKind[t]);
			else
				sprintf(cmd,"%s/MPEGGdecompress -i %s -A chr%d:%d-%d -%c",gGTLbindir,fn,chr,pos,pos+kGenomicChunk+kMaxScreenCol,tagKind[t]);
			//sleep(2);
			if (verbose)
				fprintf(stderr, "%s\n",cmd);
			f[t] = popen (cmd, "r");
			if (!f[t])
			{
				fprintf(stderr, "Error:Cannot execute command %s\n",cmd);
				return(-1);
			}
		}
	}
	else
	{
		fCnt = 1;
		sprintf(cmd,"./RLdump.pl -i %s -C %d -P %d..%d",fn,chr,pos,pos+kGenomicChunk+kMaxScreenCol);
		fprintf(stderr, "%s\n",cmd);
		//sleep(2);
		if (verbose)
			fprintf(stderr,"fetching tags\n");
		f[0] = popen (cmd, "r");
		if (!f[0])
		{
			fprintf(stderr, "Error:Cannot execute command %s\n",cmd);
			return(-1);
		}
	}
	
	// read streams
	for (t = 0; t<fCnt;t++)
	{
		while(!feof(f[t]))
		{
			int len;
			int lastP = 0;

			if (i >= gMaxTags)
			{
				fprintf(stderr, "Number of tags returned is greater than gMaxTags (%u); skipping the overflow\n", gMaxTags);
				break;
			}
			memset(tag[i].tagname,' ',kMaxTagnameLength);
			// fgets returns NULL upon encountering EOF when no characters have been read
			if (fgets(buf,kMaxLineBuf,f[t]) == NULL)
				break;
		//	fprintf(stderr,"L%6d %s\n",i,buf);
			sscanf(buf,"%s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu"
					,tag[i].tagname
					,tag[i].tagseq
					,tag[i].tagqual
					,tag[i].taginsertion
					,tag[i].taginsertionQual
					,&tag[i].tagpos
					,&tag[i].tagpair
					,&tag[i].tagLength
					,&tag[i].flag
					,&tag[i].ordinal
					);
			if (tag[i].taginsertion[0] == '=')
			{
				strcpy(tag[i].taginsertion,tag[i].tagseq);
				strcpy(tag[i].taginsertionQual,tag[i].tagqual);
			}
			len = strlen(tag[i].tagname);
			tag[i].tagname[len] = ' ';
			if (len > gMaxTagNameLength)
				gMaxTagNameLength = len;

			if (tag[i].tagpair == kIsPairedEnd)
			{
				tag[i].tagpair = i+1;
				i++;
				memset(tag[i].tagname,' ',kMaxTagnameLength);
				fgets(buf,kMaxLineBuf,f[t]);
		//		fprintf(stderr,"R%6d %s\n",i,buf);
				sscanf(buf,"%s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu"
					,tag[i].tagname
					,tag[i].tagseq
					,tag[i].tagqual
					,tag[i].taginsertion
					,tag[i].taginsertionQual
					,&tag[i].tagpos
					,&tag[i].tagpair
					,&tag[i].tagLength
					,&tag[i].flag
					,&tag[i].ordinal
					);
				if (tag[i].taginsertion[0] == '=')
				{
					strcpy(tag[i].taginsertion,tag[i].tagseq);
					strcpy(tag[i].taginsertionQual,tag[i].tagqual);
				}
				tag[i].tagpair = i-1;
				if (tag[i-1].tagpos != prev1 || tag[i].tagpos != prev2 || !(showmode & kFilterClone))
				{
					getCoverage(pos,bam,tag + i - 1,&lastP);
					getCoverage(pos,bam,tag + i,&lastP);
				}
				prev1 = tag[i-1].tagpos;
				prev2 = tag[i].tagpos;
				len = strlen(tag[i].tagname);
				tag[i].tagname[len] = ' ';
				if (len > gMaxTagNameLength)
					gMaxTagNameLength = len;
				
				unsigned int fraglen = tag[i].tagpos - tag[i-1].tagpos + tag[i].tagLength;
				if (fraglen < 1024)
					fragSizeCnt[fraglen]++;
					
			}
			else
			{
				if (tag[i].tagpos != prev1 || !(showmode & kFilterClone))
					getCoverage(pos,bam,tag + i,&lastP);
				prev1 = tag[i].tagpos;
			}
			i++;
		}
		pclose(f[t]);
		if (verbose)
			fprintf(stderr, "%d TAGS retained upon completion of %c\n",i,tagKind[t]);
	}

	bam->retained = i;
	bam->readTags = i;
	if (gMaxTagNameLength >= kMaxTagnameLength)
		gMaxTagNameLength = kMaxTagnameLength - 1;
	for (i=0; i < bam->retained; i++)
		tag[i].tagname[gMaxTagNameLength] = '\0';

	if (verbose)
		fprintf(stderr, "%d TAGS retained\n",bam->retained);
	//sleep(2);
	bam->startPos = pos;
	bam->maxCoverage = 0;
	short tem[kGenomicChunk+kMaxScreenCol];
	unsigned int tc = 0;
	for (i = 0; i < kCoverageElements * (kGenomicChunk+kMaxScreenCol); i += kCoverageElements)
	{
		unsigned int cov = bam->coverage[i] + bam->coverage[i+1] + bam->coverage[i+2] + bam->coverage[i+3] + bam->coverage[i+4] + bam->coverage[i+5] + bam->coverage[i+6];
		// put some figure.. enough to get out of noise in exome, but not too much...
		if (cov > 7)
			tem[tc++] = cov;
		//fprintf(stderr,"%u:\t%hu\t%hu\t%hu\t%hu\t%hu\t%hu\t%hu\n", pos + i / 6, bam->coverage[i], bam->coverage[i+1], bam->coverage[i+2], bam->coverage[i+3], bam->coverage[i+4], bam->coverage[i+5], bam->coverage[i+6]);
		if (cov > bam->maxCoverage)
			bam->maxCoverage = cov;
	}
	qsort(tem,tc,2,compare_short);
	bam->medianCoverage = 1;
	if (tc > 0)
		bam->medianCoverage = tem[tc/2];

	unsigned int maxCnt = 0;
	bam->mostFrequentFragmentSize = 0;
	for (i = 1; i < 1024; i++)
	{
		if (fragSizeCnt[i] > maxCnt)
		{
			maxCnt = fragSizeCnt[i];
			bam->mostFrequentFragmentSize = i;
		}
	}

	return(0);

} // getGTLchunk
// --------------------------------------------------------------------------------------------------
static int getGenomicChumk(char *fetchGenome,int chr,int pos)
{
	FILE *f = NULL;
	char cmd[kMaxFileName*4];

#ifndef ROBIN
	if (*fetchGenome != 0)
	{
#endif
		sprintf(cmd,"fetch -w 32760 %s:NC_%06d[%u..%u] | tail -1 ",fetchGenome,chr,pos+1,pos+kGenomicChunk+kMaxScreenCol); //+1 because 0 based...
		if (verbose)
			fprintf(stderr,"fetching genomic sequence\n");
#ifndef ROBIN
	}
	else
	{
		if (gGTLbindir[0] == 0)
			sprintf(cmd,"GTLfetch -g %s -C %d -P %u..%u",gGTLgenome,chr,pos+1,pos+kGenomicChunk+kMaxScreenCol); //+1 because 0 based...
		else
			sprintf(cmd,"%s/GTLfetch -g %s -C %d -P %u..%u",gGTLbindir,gGTLgenome,chr,pos+1,pos+kGenomicChunk+kMaxScreenCol); //+1 because 0 based...
		if (verbose)
			fprintf(stderr,"fetching genomic sequence using GTLfetch\n%s\n",cmd);
	}
#endif
	f = popen (cmd, "r");
	if (!f)
	{
		fprintf(stderr, "Error:Cannot execute command %s\n",cmd);
		return(1);
	}
	fgets(genome,kMaxLineBuf,f);
	pclose(f);
	if (gShowAnnotation != NULL)
	{
		sprintf(cmd,"%s/P_NC_%06d.txt",gShowAnnotation,chr);
		if (verbose)
			fprintf(stderr,"fetching plus strand annotation '%s'\n",cmd);
		f = fopen (cmd, "r");
		if (!f)
		{
			fprintf(stderr, "Error:Cannot open %s\n",cmd);
			return(1);
		}
		fseek(f,pos,SEEK_SET);
		fgets(gAnnotationFwd,kMaxLineBuf,f);
		fclose(f);
		sprintf(cmd,"%s/M_NC_%06d.txt",gShowAnnotation,chr);
		if (verbose)
			fprintf(stderr,"fetching minus strand annotation '%s'\n",cmd);
		f = fopen (cmd, "r");
		if (!f)
		{
			fprintf(stderr, "Error:Cannot open %s\n",cmd);
			return(1);
		}
		fseek(f,pos,SEEK_SET);
		fgets(gAnnotationRev,kMaxLineBuf,f);
		fclose(f);
		sprintf(cmd,"%s/N_NC_%06d.txt",gShowAnnotation,chr);
		if (verbose)
			fprintf(stderr,"fetching gene name annotation '%s'\n",cmd);
		f = fopen (cmd, "r");
		if (!f)
		{
			fprintf(stderr, "Error:Cannot open %s\n",cmd);
			return(1);
		}
		fseek(f,pos,SEEK_SET);
		fgets(gAnnotationGene,kMaxLineBuf,f);
		fclose(f);
#if 0
		unsigned int i=0;
		while(gAnnotationFwd[i] != 0)
		{
			if (gAnnotationFwd[i] == '}')
				gAnnotationFwd[i] = '>';
			else if (gAnnotationFwd[i] == '.')
				gAnnotationFwd[i] = ' ';
			i++;
		}
		i=0;
		while(gAnnotationRev[i] != 0)
		{
			if (gAnnotationRev[i] == '{')
				gAnnotationRev[i] = '<';
			else if (gAnnotationRev[i] == '.')
				gAnnotationRev[i] = ' ';
			i++;
		}
		i=1;
		while(gAnnotationGene[i] != 0)
		{
			if ((gAnnotationGene[i-1] == '-') && (gAnnotationGene[i] == '-'))
			{
				gAnnotationGene[i-1] = ' ';
				gAnnotationGene[i] = ' ';
			}
			else if ((gAnnotationGene[i-1] == '.') && (gAnnotationGene[i] == '.'))
			{
				gAnnotationGene[i-1] = ' ';
				gAnnotationGene[i] = ' ';
			}
			i++;
		}
#endif
	}

	return(0);

} /* getGenomicChumk */
// --------------------------------------------------------------------------------------------------
static int scrollView( int *availFromPos,unsigned int *curpos, int amount,int row)
{
	*curpos += amount;
	if ((*curpos < *availFromPos) || (*curpos > (*availFromPos+kGenomicChunk-kMaxScreenCol)))
	{
		*availFromPos = (*curpos > kHalfGenomicChunk) ? *curpos - kHalfGenomicChunk : 0;
		printf("\033[%d;1f\033[JFetching genomic chunk starting from %d",row,*availFromPos);
		fflush(stdout);
		return(1);
	}
	return(0);

} /* scrollView */
// --------------------------------------------------------------------------------------------------
static void DisplayPatientID(BAMFILE *bam)
{
	int from = strlen(bam->patientname)-16;
	if (from < 0)
		from = 0;

	strncpy(bam->displayname,&bam->patientname[from],12);


} /* DisplayPatientID */
// --------------------------------------------------------------------------------------------------
static void DisplayBAMname(BAMFILE *bam)
{
	int from = strlen(bam->fn)-16;
	if (from < 0)
		from = 0;

	strncpy(bam->displayname,&bam->fn[from],12);

} /* DisplayBAMname */
// --------------------------------------------------------------------------------------------------

static int GetPatientsHighlight(char *fn,PATIENT *p)
{
	FILE *f = NULL;
	char cmd[kMaxFileName];
	char linbuf[kMaxLineBuf];
	int cnt = 0;



	f = fopen (fn, "r");
	if (!f)
	{
		fprintf(stderr, "Error:Cannot execute command %s\n",cmd);
		return 0;
	}

	fgets(linbuf,kMaxLineBuf,f); // skip header.
	while(!feof(f))
	{
		fgets(linbuf,kMaxLineBuf,f);
		if (feof(f))
			break;
		sscanf(&linbuf[0],"%s\t%s\t%d",p[cnt].bamfile,p[cnt].patient,&p[cnt].color);
		if (++cnt == (kMaxPatients-1))
			break;
#if 0
		int pid = 0;
		while(p[pid].bamfile[0] != 0)
		{
			if (strcmp(patient,p[pid].patient) == 0)
			{
				p[pid].color = color;
				break;
			}
			pid++;
		}
#endif
	}
	fclose(f);

	return(cnt);

} /* GetPatientsHighlight */
// --------------------------------------------------------------------------------------------------

static void
loadSNP(NAV *n, FILE *f, const char *snpfile)
{
	char linbuf[kMaxLineBuf];

	while(fgets(linbuf,kMaxLineBuf,f) != NULL)
	{
		char tmpchr[32];
		int pos,bg,fg;
		int cnt = sscanf(&linbuf[0],"%s\t%u\t%u\t%u",tmpchr,&pos,&bg,&fg);
		if (cnt < 2)
			continue;
		if (n->nb == n->max)
		{
			n->max += 1024;
			n->loc = realloc(n->loc, n->max * sizeof(SNPTABLE));
			if (n->loc == NULL)
			{
				fprintf(stderr, "Error: realloc failed while reading %s\n",snpfile);
				fclose(f);
				n->nb = 0;
				n->max = 0;
				return;
			}
		}
		SNPTABLE *cs = n->loc + n->nb;
		n->nb += 1;
		cs->pos = pos;
		cs->fgCode = (cnt < 4) ? 37 : fg; // default to white
		cs->bgCode = (cnt < 3) ? 44 : bg; // default to blue
		// FIXME - should get chr name from the .cfg file if available
		if ((tmpchr[0] == 'c') && (tmpchr[1] == 'h') && (tmpchr[2] == 'r'))
		{
			if (tmpchr[3] == 'X') cs->chr = 23;
			else if (tmpchr[3] == 'Y') cs->chr = 24;
			else cs->chr = atoi(&tmpchr[3]);
		}
		else
		{
			if (tmpchr[0] == 'X') cs->chr = 23;
			else if (tmpchr[0] == 'Y') cs->chr = 24;
			else cs->chr = atoi(tmpchr);
		}
	}
}
// --------------------------------------------------------------------------------------------------
static int
loadSNPfile(NAV *n, char *snpfile)
{
	if (verbose)
		fprintf(stderr, "Reading %s\n",snpfile);

	FILE *f = fopen (snpfile, "r");
	if (!f)
	{
		fprintf(stderr, "Error:Cannot open %s\n",snpfile);
		return 1;
	}
	loadSNP(n, f, snpfile);
	fclose(f);
	return 0;
}
// --------------------------------------------------------------------------------------------------
static int readGTL(BAMFILE *bam,TAG *tag,int chr,int pos,int realign)
{
	int err;

	//fprintf(stderr, "readGTL: %s format %d chr%d:%d\n",bam->fn,bam->fileformat,chr,pos);
	//sleep(1);

	if (pos > 1000)
		pos -= 1000;
	else
		pos = 0;

	bam->retained = 0;
	bam->readTags = 0;
	err = getGTLchunk(bam->fn,chr,pos,bam,tag,realign);
	if (err < 0)
	{
		fprintf(stderr, "[ADVIEW] getGTLchunk error\n");
		return(0);
	}

	return(chr);

} // readGTL
// --------------------------------------------------------------------------------------------------
static int readBAM(BAMFILE *bam,TAG *tag,int chr,int pos)
{
	int curChr;
	char chrToProcess[32];


	//fprintf(stderr, "FILENAME: %s format %d (GTL is %d)\n",bam->fn,bam->fileformat,kIsGTL);
	//sleep(1);
#ifndef ROBIN
	if (bam->fileformat == kIsGTL || bam->fileformat == kIsTB || bam->fileformat == kIsMPEGG)
	{
		return(readGTL(bam,tag,chr,pos,0));
	}
#endif


	// in fact instead of 1000 we should probably start at pos - kMaxTagLength to make sure RNAseq (with introns) are fully read.
	// this was an arbitrary value roughly equal to tagSize (kind if 128nt or so) + max screen width would always be true...

#if 0
	if (chr <= 22)
		sprintf(chrToProcess,"%s%d:%d",gNoChrInBAM ? "" : "chr",chr,(pos > 1000) ? pos-1000 : 0);
	else if (chr == 23)
		sprintf(chrToProcess,"%sX:%d",gNoChrInBAM ? "" : "chr",(pos > 1000) ? pos-1000 : 0);
	else if (chr == 24)
		sprintf(chrToProcess,"%sY:%d",gNoChrInBAM ? "" : "chr",(pos > 1000) ? pos-1000 : 0);
#endif
	sprintf(chrToProcess,"%s:%d",(gNoChrInBAM && chrPos[chr].name[0] == 'c' && chrPos[chr].name[1] == 'h' && chrPos[chr].name[2] == 'r') ? chrPos[chr].name + 3 : chrPos[chr].name,(pos > 1000) ? pos-1000 : 0);


	/* --------- read sequencing data from BAM file ----------------*/

	if (verbose)
	{	fprintf(stderr, "processing bam\n");fflush(stderr); }
	bam->retained = 0;
	bam->readTags = 0;
	bam->iter = sam_itr_querys(bam->idx, bam->header, chrToProcess); // parse a region in the format like `chr2:100-200'
	if (bam->iter == NULL)
	{ // reference name is not found
		fprintf(stderr, "region \"%s\" specifies an unknown reference name.\n", chrToProcess);fflush(stderr);
		return(0);
	}


	bam->b = bam_init1();
	curChr = sam_readADNI(bam,bam->b,tag,pos,pos + kGenomicChunk + 1000);  // not sure we absolutely need the +1000 here
	bam_destroy1(bam->b);
	if (curChr < -1)
	{
		fprintf(stderr, "[ADVIEW] truncated file.\n");
		return(0);
	}

	return(curChr);

} /* readBAM */
// --------------------------------------------------------------------------------------------------
static void *readBAMthread(void *ep)
{
	((EXECUTIONPLAN*)ep)->rslt = readBAM(((EXECUTIONPLAN*)ep)->bam,((EXECUTIONPLAN*)ep)->tag,((EXECUTIONPLAN*)ep)->chr,((EXECUTIONPLAN*)ep)->pos);
	pthread_exit(NULL);
}
// --------------------------------------------------------------------------------------------------
static int changeActiveScreenBam(BAMFILE *bam, TAG *tag, int newbam, char* bamfiles, char *dir, int chr, int pos,PATIENT *patient)
{
	int err;

	closeBAM(bam);
	bam->curbam = newbam;
	sprintf(bam->fn,"%s/%s.bam",dir,&bamfiles[bam->curbam*kMaxBamName]);
	if (access(bam->fn,R_OK) != 0)
		sprintf(bam->fn,"%s/%s",dir,&bamfiles[bam->curbam*kMaxBamName]);
	strcpy(bam->ADNIbamname,&bamfiles[bam->curbam*kMaxBamName]);
	err = openBAM(bam,patient);
	if (err)
		return(-1);
	if (readBAM(bam,tag,chr,pos) < 0)
		return(-1);

	return(0);

} /* changeActiveScreenBam */
// --------------------------------------------------------------------------------------------------
int getNextPatientKind(int *list,int cur,int howManyTimes)
{
	int k;
	int hit = -1;

		int i = 0;
		while (list[i] != -1)
		{

			if (list[i] > cur)
			{
				break;
			}
			i++;
		}
		if (list[i] == -1)
			i = 0;

		for (k = 0; k < (howManyTimes-1); k++)
		{
			i++;
			if (list[i] == -1)
				i = 0;
		}

		hit = list[i];
		return(hit);

}	/* getNextPatientKind */
// --------------------------------------------------------------------------------------------------
int getPreviousPatientKind(int *list,int cur,int howManyTimes)
{
	int k;
	int hit = -1;
	int lastPatient = 0;

		while (list[lastPatient] != -1)
		{
			lastPatient++;
		}
		lastPatient--;
		if (lastPatient < 0)
			return(0);

		int i = lastPatient;
		while (i != -1)
		{

			if (list[i] < cur)
			{
				break;
			}
			i--;
		}
		if (i == -1)
			i = lastPatient;

		for (k = 0; k < (howManyTimes-1); k++)
		{
			i--;
			if (i == -1)
				i = lastPatient;
		}

		hit = list[i];
		return(hit);

}	/* getPreviousPatientKind */
// --------------------------------------------------------------------------------------------------

static int selectbam(char *bamfiles,int bamcnt,BAMFILE *bam,int activeScreen,int nrows, int ncol,int colsize,PATIENT *p)
{
	int x,i;
	int shown = 0;
	int from = 0;
	char selection[64];
	int pos = 0;
	int c;
	int k;
	int choice = bam[activeScreen].curbam;


	char initBGcolor[16];
	char endBGcolor[8];

	redraw:

		printf("\033[H\033[J");
		shown = 0;
		for (x =0; x < ncol; x++)
		{
			for (i =0; i< nrows; i++)
			{
				int idx = (x*nrows+i+from);
				if (idx >= bamcnt)
					goto last;
				printf("\033[%d;%df",i+1,x*colsize+1);

				int color = 0;
				int pid = 0;
				char pn[64];
				while(p[pid].bamfile[0] != 0)
				{
					if (strcmp(&bamfiles[idx*kMaxBamName],p[pid].bamfile) == 0)
					{
						color = p[pid].color;
						strcpy(pn,p[pid].patient);
						break;
					}
					pid++;
				}

				/* check color */
				switch(color)
				{
					case 0: initBGcolor[0] = 0;  endBGcolor[0] = 0;  break;
					case 1: strcpy(initBGcolor,"\033[1;42m\033[1;37m");  strcpy(endBGcolor,"\033[0m");  break;  // green
					case 2: strcpy(initBGcolor,"\033[1;43m\033[1;37m");  strcpy(endBGcolor,"\033[0m");  break;  // yellow
					case 3: strcpy(initBGcolor,"\033[1;41m\033[1;37m");  strcpy(endBGcolor,"\033[0m");  break;  // red
				}


				// check if visible in an otherscreen
				for (k=0; k<gScreenCnt;k++)
				{
					if ((shown+from) == bam[k].curbam)
					{
						if (bam[k].curbam == bam[activeScreen].curbam)
							printf("%s%4u)%s \033[1;31m%s\033[0m\n",initBGcolor,(shown+from+1),endBGcolor,(showmode & kDisplayPatientID) ? pn : &bamfiles[idx*kMaxBamName]); // red
						else
							printf("%s%4u)%s \033[1;33m%s\033[0m\n",initBGcolor,(shown+from+1),endBGcolor,(showmode & kDisplayPatientID) ? pn : &bamfiles[idx*kMaxBamName]); // orange
						break;
					}
				}
				if (k == gScreenCnt)
				{
					printf("%s%4u)%s %s\n",initBGcolor,(shown+from+1),endBGcolor,(showmode & kDisplayPatientID) ? pn : &bamfiles[idx*kMaxBamName]); // normal
				}
				shown++;
			}
		}

	last:;

		while(1)
		{
			read(0,&c,1);
			switch((char)c)
			{
				case 'n':
					if (from < (bamcnt -1))
						from += (nrows*ncol);
					goto redraw;

				case 'p':
					if (from > 0)
						from -= (nrows*ncol);
					goto redraw;

				case 'd':
					if (showmode & kDisplayPatientID)
					{
						showmode &= ~kDisplayPatientID;
					}
					else
					{
						showmode |= kDisplayPatientID;
					}
					goto redraw;

				case '0':
				case '1':
				case '2':
				case '3':
				case '4':
				case '5':
				case '6':
				case '7':
				case '8':
				case '9':
					selection[pos++] = (char)c;
					putchar((char)c);
					fflush(stdout);
					break;

				case 0x08:
					if (pos > 0)
						pos--;
					selection[pos] = 0;
					break;

				case '\n':
					printf("\033[J");
					selection[pos] = 0;
					choice = atoi(selection) - 1;  // user is 1 based, but not us !
					if ((choice < 0) || (choice >= bamcnt))
					{
						printf("Invalid choice. Will do nothing.");
						fflush(stdout);
						choice = bam[activeScreen].curbam;
						sleep(1);
					}
					goto done;

				case 0x1b:
				case 'q':
					printf("\033[J");
					goto done;
			}
		}

	done:;
		return(choice);

} /* selectbam */
// --------------------------------------------------------------------------------------------------
static void SIGQUITcatch (int sig)
{
	printf("\033[JPlease press Q to quit cleanly.");
	fflush(stdout);
}
// --------------------------------------------------------------------------------------------------

//chr	len	virtChr	offset	AC	SAMname
//1	249250621	1	0	1	chr1
//2	243199373	2	0	2	chr2
//3	198022430	3	0	3	chr3
//16	90354753	0	0	16	chr16
//17	81195210	0	91000000	17	chr17
//27	171823	1	251000000	NC_007605	chrNC_007605
static void
load_chr_pos(const char *name)
{
	size_t l = strlen(name);
	char fname[l + 1];
	strcpy(fname,name);
	strcpy(fname + l - 4, ".cfg");
	FILE *f = fopen(fname,"r");
	char s[512];
	char n[32];
	while (fgets(s,512,f) != NULL)
	{
		unsigned int chr,len,vc,offset;
		if (sscanf(s,"%u %u %u %u %*s %s",&chr,&len,&vc,&offset,n) == 5)
		{
			nbChrPos = chr + 1;
			chrPos[chr].len = len;
			chrPos[chr].offset = (vc << 28) + offset;
			strcpy(chrPos[chr].name,n);
			if (verbose)
				fprintf(stderr,"got %u %u %u %s\n",chr,chrPos[chr].len,chrPos[chr].offset,chrPos[chr].name);
		}
	}
	fclose(f);
}
// --------------------------------------------------------------------------------------------------
static int
find_closest_match(NAV *n, int chr, int pos)
{
	int best = -1;
	int i;
	for (i = 0; i < n->nb; i++) {
		if (n->loc[i].chr == chr) {
			if (n->loc[i].pos == pos) {
				return i;
			} else {
				if (best == -1)
					best = i;
				if (abs(n->loc[i].pos - pos) < abs(n->loc[best].pos - pos))
					best = i;
			}
		}
	}
	return best;
}
// --------------------------------------------------------------------------------------------------

int main_ADinteractive(int argc, char *argv[])
{
	DIR *dirlist = NULL;
	char *bamfiles = NULL;
	TAG *tag[kMaxVirutalScreen];
	BAMFILE bam[kMaxVirutalScreen];
	VIRTUALSCREEN screenpanel[kMaxVirutalScreen];
	char dir[kMaxFileName];
	char snpfile[kMaxFileName];
	char recsnpfile[kMaxFileName];
	char highlightfile[kMaxFileName];
	char fetchGenome[512];
	PATIENT patient[kMaxPatients];
	int patientcnt = 0;
	int listCN[kMaxDir+1];
	int listAD[kMaxDir+1];
	int listMCI[kMaxDir+1];

	int bamcnt = 0;
	int activeScreen = 0;
	int desiredScreenCount = kMaxVirutalScreen;

	unsigned int linkedSNP1;
	unsigned int linkedSNP2;
	char linkedSNPfn[kMaxFileName];

	int chr = -1;
	int availFromPos = 0;

	unsigned int curpos = 0;
	struct winsize w;
	struct termios old,new;
	unsigned int line = 1;
	int maxFileNameWidth = 0;
	int updateGenomeSeq = 1;
	int midscreen;
	int i;
	int c;
	char *preExecute = NULL;
	int err = 0;
	alignerConfigFile[0] = '\0';


#ifdef DEBUGOUTPUT
debug  = fopen("/tmp/dbg.txt","w");
#endif

//	setvbuf(stdout, buf, _IOLBF , 8192);

	/* ------  initialize a few variables ----- */

	for (i = 0; i < kMaxVirutalScreen; i++)
	{
		tag[i] = NULL;
		bam[i].in = 0;
		bam[i].header = NULL;
		bam[i].idx = NULL;
		bam[i].iter = NULL;
		bam[i].curbam = -1;
		bam[i].patientname[0] = 0;
		bam[i].ADNIbamname[0] = 0;
		bam[i].displayname[0] = 0;
		bam[i].info[0] = 0;
		bam[i].fn[0] = 0;
		bam[i].color = 0;
		bam[i].fileformat = 0;
	}
	dir[0] = 0;
	snpfile[0] = 0;
	recsnpfile[0] = 0;
	highlightfile[0] = 0;
	opterr = 0;
#ifdef ROBIN
	strcpy(fetchGenome,"hg19");
#else
	fetchGenome[0] = 0;
#endif
	listAD[0] = -1;
	listMCI[0] = -1;
	listCN[0] = -1;
	linkedSNPfn[0] = 0;
	{
		char *s = getenv("GTLBINDIR");
		if (s != NULL)
		{
			unsigned int i;
			strncpy(gGTLbindir,s,511);
			s[511] = 0;
			for (i = 0; i < 512; i++)
			{
				if (s[i] != '/' && s[i] != '.' && s[i] != '-' && s[i] != '_' && s[i] != '~' && !isalnum(s[i]))
				{
					s[i] = 0;
					break;
				}
			}
		}
	}
	/* --------- process arguments */

	signal (SIGINT, SIGQUITcatch);
	tcgetattr(0,&old);

#ifdef ROBIN
	gMaxTags = kMaxTagsXScreens ;
	while ((c = getopt (argc, argv, "1:c:p:SD:qm:l:")) != -1)
#else
	while ((c = getopt (argc, argv, "C1:2:3:4:5:6:7:8:9:a:d:n:v:c:p:r:s:g:G:h:A:B:L:SD:zU:e:")) != -1)
#endif
	switch (c)
	{
		case 'C':				gNoCompare = 1; break;
		case 'z':				gNoChrInBAM = 1; break;
		case 'A': 			sscanf(optarg,"%d",&linkedSNP1);  linkedSNP1--; break; // convert zero based
		case 'B': 			sscanf(optarg,"%d",&linkedSNP2);  linkedSNP2--; break;
		case 'L': 			strcpy(linkedSNPfn,optarg); break;

		case 'a':
			gShowAnnotation = optarg;
		break;

		case 'e':
			preExecute = optarg;
		break;

		case '1':
			strcpy(bam[0].fn,optarg);
			if (gNbVirtualScreen < 1)
			gNbVirtualScreen = 1;
		break;

		case '2':
			strcpy(bam[1].fn,optarg);
			if (gNbVirtualScreen < 2)
			gNbVirtualScreen = 2;
		break;

		case '3':
			strcpy(bam[2].fn,optarg);
			if (gNbVirtualScreen < 3)
			gNbVirtualScreen = 3;
		break;

		case '4':
			strcpy(bam[3].fn,optarg);
			if (gNbVirtualScreen < 4)
			gNbVirtualScreen = 4;
		break;

		case '5':
			strcpy(bam[4].fn,optarg);
			if (gNbVirtualScreen < 5)
			gNbVirtualScreen = 5;
		break;

		case '6':
			strcpy(bam[5].fn,optarg);
			if (gNbVirtualScreen < 6)
			gNbVirtualScreen = 6;
		break;

		case '7':
			strcpy(bam[6].fn,optarg);
			if (gNbVirtualScreen < 7)
			gNbVirtualScreen = 7;
		break;

		case '8':
			strcpy(bam[7].fn,optarg);
			if (gNbVirtualScreen < 8)
			gNbVirtualScreen = 8;
		break;

		case '9':
			strcpy(bam[8].fn,optarg);
			if (gNbVirtualScreen < 9)
			gNbVirtualScreen = 9;
		break;

		case 'd':
			strcpy(dir,optarg);
		break;

		case 'g':
			strcpy(fetchGenome,optarg);
			break;

		case 'U':
			strcpy(alignerConfigFile,optarg);
		break;

		case 'G':
			strcpy(gGTLgenome,optarg);
			load_chr_pos(gGTLgenome);
		break;

		case 'r':
			strcpy(recsnpfile,optarg);
		break;

		case 's':
			strcpy(snpfile,optarg);
		break;

		case 'h':
			strcpy(highlightfile,optarg);
		break;

		case 'c':
			sscanf(optarg,"%d",&chr);
		break;

		case 'p':
			sscanf(optarg,"%d",&curpos);
			availFromPos = (curpos > kHalfGenomicChunk) ? curpos - kHalfGenomicChunk : 0;
		break;

#ifdef ROBIN
		case 'l':
			sscanf(optarg,"%d",&gRobinPrintEveryNthLine);
		break;

		case 'q':
			showmode |= kShowQuality;
		break;

		case 'm':
			sscanf(optarg,"%u",&gMaxTags);
		break;
#endif

		case 'n':
			sscanf(optarg,"%d",&desiredScreenCount);
			if (desiredScreenCount > kMaxVirutalScreen)
				desiredScreenCount = kMaxVirutalScreen;
			if (desiredScreenCount < 1)
				desiredScreenCount = 1;
			if (gNbVirtualScreen < desiredScreenCount)
				gNbVirtualScreen = desiredScreenCount;
		break;

		case 'v':
			sscanf(optarg,"%d",&verbose);
		break;

		case 'S':
			singleEnd = 1;
		break;

		case 'D':
			sscanf(optarg,"%d",&directDump);
			if (directDump > (kMaxScreenCol-1))
				directDump = (kMaxScreenCol-1);
			if (directDump < 1)
				directDump = 1;
		break;
	}


	if ( ((bam[0].fn[0] == 0) && (dir[0] == 0)) || (chr == -1) || (curpos == 0) )
	{
		usage();
		err = 1;
		goto bail;
	}

	/* ------  allocate memory to read tags ----- */

	if (verbose)
	{	fprintf(stderr, "Allocating memory\n");fflush(stderr); }


#ifdef ROBIN
	err = posix_memalign((void **)&tag[0], 64, gMaxTags*sizeof(TAG));
	if (err != 0)
	{
		fprintf(stderr, "Error: 1 cannot allocate memory\n");
		err = 1;
		goto bail;
	}
#else
	gMaxTags = kMaxTagsXScreens / gNbVirtualScreen;
	for (i = 0; i < gNbVirtualScreen; i++)
	{
		err = posix_memalign((void **)&tag[i], 64, gMaxTags*sizeof(TAG));
		err += posix_memalign((void **)&screenpanel[i].virtualscreen, 64, (12*kMaxScreenCol*kVirtScreenMaxLines)* sizeof(char)); // some space for color chars.
		if (err != 0)
		{
			fprintf(stderr, "Error: 1 cannot allocate memory\n");
			err = 1;
			goto bail;
		}
	}
#endif

	/* ------  allocate memory to read directory ----- */

	err = posix_memalign((void **)&bamfiles, 64, kMaxDir*kMaxBamName*sizeof(char));
	if (err != 0)
	{
		fprintf(stderr, "Error: 2 cannot allocate memory\n");
		err = 1;
		goto bail;
	}

	/* --------- get patient info ----------------*/


	if (highlightfile[0])
		patientcnt = GetPatientsHighlight(highlightfile,patient);

#ifndef ROBIN
	/* -- read directory -- */
	if (dir[0] != 0)
	{
		if (verbose)
			fprintf(stderr, "Finding indexed bam files in %s\n",snpfile);
		dirlist = opendir(dir);
		if (dirlist)
		{
			struct dirent *tmp;

			while((tmp=readdir(dirlist))!=NULL)
			{
				int l = strlen(tmp->d_name);
				if (l > maxFileNameWidth)
					maxFileNameWidth = l;
				if (l > 8)
				{
					if (tmp->d_name[l-4] == '.' && tmp->d_name[l-3] == 'b' && tmp->d_name[l-2] == 'a' && tmp->d_name[l-1] == 'i')
					{
						memcpy(&bamfiles[bamcnt*kMaxBamName],tmp->d_name,l-4);
						bamfiles[bamcnt*kMaxBamName + l - 8] = '\0'; // strip .bam.bai
						bamcnt++;
					}
				}
				else if (tmp->d_type == DT_DIR || tmp->d_type == DT_LNK)
				{
					char tfn[1024];
					sprintf(tfn,"%s/%s/chr1_i.txt",dir,tmp->d_name);
					if (access(tfn,R_OK) == 0)
					{
						strcpy(&bamfiles[bamcnt*kMaxBamName],tmp->d_name);
						bamcnt++;
					}
				}
				if (bamcnt == kMaxDir)
					break;
			}
			closedir(dirlist);
			sortbamNames(bamfiles,bamcnt);
			if (patientcnt > 0)
				getPatientKindLists(bamfiles,bamcnt,patient,listCN,listMCI,listAD);

		}
		if (bam[0].fn[0] == 0)
		{
			if (bamcnt == 0)
			{
				fprintf(stderr, "Error: no indexed bam found in %s/\n",dir);
				err = 1;
				goto bail;
			}
			// open up to kMaxVirutalScreen files
			for (i=0; i<desiredScreenCount;i++)
			{
				if (i < bamcnt)
				{
					sprintf(bam[i].fn,"%s/%s.bam",dir,&bamfiles[i*kMaxBamName]);
					if (access(bam[i].fn,R_OK) != 0)
						sprintf(bam[i].fn,"%s/%s",dir,&bamfiles[i*kMaxBamName]);
					strcpy(bam[i].ADNIbamname,&bamfiles[i*kMaxBamName]);
					bam[i].curbam = i;
				}
			}
		}
	}

	/* --------- get snp file ----------------*/

	if (snpfile[0] != 0)
	{
		err = loadSNPfile(&SNPpos, snpfile);
		if (err)
			goto bail;
		curNPpos = &SNPpos; // so that n and p cycle through the SNPs
	}

	/* --------- get recorded snp file ----------------*/

	if (recsnpfile[0] != 0)
	{
		err = loadSNPfile(&recSNPpos, snpfile);
		if (err)
			goto bail;
	}
#endif
	/* --------- open bam(s) file ----- */

	err = openBAM(&bam[0],patient);
	if (err)
		goto bail;

	gScreenCnt = 1;
	for (i = 1; i < gNbVirtualScreen; i++)
	{
		if (bam[i].fn[0])
		{
			showmode |= kSplitScreen;
			gScreenCnt++;
			err = openBAM(&bam[i],patient);
			if (err)
				goto bail;
		}
	}
	// clear screen.
	if (!directDump)
		printf("\033[H\033[J\n");

// --------------------------------------------------------------------------------------------------
//                MAIN EVENT loop
// --------------------------------------------------------------------------------------------------


	/* --------- move to chr position ----------------*/

changechr:;

	/* --------- get genomic sequence ----------------*/

	if (updateGenomeSeq)
	{
#ifdef DEBUGOUTPUT
fprintf(debug,"changechr chr=%u; availFromPos=%u\n",chr,availFromPos);
#endif
		err = getGenomicChumk(fetchGenome,chr,availFromPos);
		if (err)
			goto bail;
	}


	/* --------- read sequencing data from BAM file ----------------*/

	if (showmode & kSplitScreen)
	{
		EXECUTIONPLAN ep[gNbVirtualScreen];
		pthread_t thread[gNbVirtualScreen];

		for (i = 0; i < gScreenCnt; i++)
		{
			ep[i].bam = &bam[i];
			ep[i].tag = tag[i];
			ep[i].chr = chr;
			ep[i].pos = availFromPos;
			ep[i].screenpanel = i;
			if (pthread_create (&thread[i], NULL, &readBAMthread, &ep[i]) == 0) {
				if (verbose > 1) { fprintf(stderr,"Launching thread %u\n",i) ;}
			}
		}

		for (i = 0; i < gScreenCnt; i++)
			if (pthread_join (thread[i], NULL) == 0) { if (verbose > 1) { fprintf(stderr,"Finished thread %u \n",i) ;} }

		for (i = 0; i < gScreenCnt; i++)
		{
			if (ep[i].rslt < 0)
				goto bail;
		}

		compareAlignments(tag[0],tag[1],bam[0].retained,bam[1].retained);

	}
	else
	{
		err = readBAM(&bam[0],tag[0],chr,availFromPos);
		if (err < 0)
		{
			fprintf(stderr,"EXIT:chr%d\n",err);
			goto bail;
		}
	}


	if (directDump)
	{
		// for Robin and Lou for HUG web browser view of amplicons
		showmode |= kShowTagName;
		printchunkForRobin(tag[0],bam[0].retained,chr,availFromPos,availFromPos+kHalfGenomicChunk,directDump);
		exit(0);
	}

	if (linkedSNPfn[0] != 0)
	{
			showmode |= kMaskLowQuality;
			goto autoprocess;
	}
	/* ------- set terminal display accordingly to our needs ------*/

#if defined __APPLE__
// whatever - just define something already in the list
#define IUCLC IXANY
#define OLCUC OCRNL
#define XCASE ECHONL
#endif
	// CI - grmp ... ???  Apparently, the terminal processing gets garbled on some systems, so try to force something sensible
	old.c_iflag &= ~(IGNBRK | PARMRK | INPCK | ISTRIP | INLCR | IGNCR | IUCLC | IXANY | IXOFF);
	old.c_iflag |= (BRKINT | IGNPAR | ICRNL | IXON | IMAXBEL);
	old.c_oflag &= ~(OLCUC | OCRNL | ONOCR | ONLRET | OFILL | OFDEL | NLDLY | CRDLY | TABDLY | BSDLY | VTDLY | FFDLY);
	old.c_oflag |= (OPOST | ONLCR);
	old.c_cflag &= ~(CSIZE | CSTOPB | PARENB | PARODD | CLOCAL);
	old.c_cflag |= (CS8 | CREAD);
	old.c_lflag &= ~(ECHONL | NOFLSH | TOSTOP | ECHOPRT | XCASE);
	old.c_lflag |= (ISIG | ICANON | ECHO | ECHOE | ECHOK | ECHOCTL | ECHOKE);

	new = old;
	new.c_lflag &= ~ICANON;
	new.c_lflag &= ~ECHO;
	new.c_cc[VMIN] = 1;
	new.c_cc[VTIME] = 0;
	tcsetattr(0,TCSANOW,&new);

	while (keeprunning)
	{
		ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);

		// print alignment
		midscreen = (w.ws_col-kScreenPanelMargin) / 2;

		if (showmode & kZoomActivePanel) // fool the system, printing always the active panel in the first one, as a full screen :-)
		{
			printchunk(tag[activeScreen],bam + activeScreen,chr,line,availFromPos,curpos,(w.ws_col-kScreenPanelMargin*1)/1,w.ws_row-1,&screenpanel[0],bam[activeScreen].color);
			flushscreen(screenpanel,1);
		}
		else
		{
			midscreen /= gScreenCnt;
			for (i = 0; i < gScreenCnt; i++)
				printchunk(tag[i],bam + i,chr,line,availFromPos,curpos,(w.ws_col-kScreenPanelMargin*gScreenCnt)/gScreenCnt,w.ws_row-1,&screenpanel[i],bam[i].color);
			flushscreen(screenpanel,gScreenCnt);
		}



		// highlight all known SNPs
		if (showmode & kShowSNP)
		{
			int k;
			unsigned int voff = 3;
			if (curpos - bam[0].startPos <= kGenomicChunk && (showmode & kShowCoverage) && !(showmode & kZoomCoverage))
				voff += ((w.ws_row - 1) > 34) ? 9 : 3;

			if (showmode & kZoomActivePanel) // fool the system, printing always the active panel in the first one, as a full screen :-)
			{
				for (k=0; k< SNPpos.nb; k++)
				{
					SNPTABLE *cs = SNPpos.loc + k;
					if ((cs->chr == chr) && (cs->pos > curpos) && (cs->pos < curpos+(w.ws_col-kScreenPanelMargin*1)/1))
					{
						printf("\033[%u;%df\033[1;%um\033[1;%um%c\033[0m\n",
							voff,
							(cs->pos-curpos)+0*((w.ws_col-kScreenPanelMargin*1)/1+kColOffset)+kColOffset,
							cs->bgCode,
							cs->fgCode,
							genome[cs->pos-availFromPos-1]); // +kColOffset because of offset line mnumbers +1  -1 because genome is zero base ?
					}
				}
				for (k=0; k< recSNPpos.nb; k++)
				{
					SNPTABLE *cs = recSNPpos.loc + k;
					if ((cs->chr == chr) && (cs->pos > curpos) && (cs->pos < curpos+(w.ws_col-kScreenPanelMargin*1)/1))
					{
						printf("\033[%u;%df\033[1;%um\033[1;%um%c\033[0m\n",
							voff,
							(cs->pos-curpos)+0*((w.ws_col-kScreenPanelMargin*1)/1+kColOffset)+kColOffset,
							cs->bgCode,
							cs->fgCode,
							genome[cs->pos-availFromPos-1]); // +kColOffset because of offset line mnumbers +1  -1 because genome is zero base ?
					}
				}
			}
			else
			{
				for (k=0; k< SNPpos.nb; k++)
				{
					SNPTABLE *cs = SNPpos.loc + k;
					if ((cs->chr == chr) && (cs->pos > curpos) && (cs->pos < curpos+(w.ws_col-kScreenPanelMargin*gScreenCnt)/gScreenCnt))
					{
						for (i = 0; i < gScreenCnt; i++)
							printf("\033[%u;%df\033[1;%um\033[1;%um%c\033[0m\n",
								voff,
								(cs->pos-curpos)+i*((w.ws_col-kScreenPanelMargin*gScreenCnt)/gScreenCnt+kColOffset)+kColOffset,
								cs->bgCode,
								cs->fgCode,
								genome[cs->pos-availFromPos-1]); // +kColOffset because of offset line mnumbers +1  -1 because genome is zero base ?
					}
				}
				for (k=0; k< recSNPpos.nb; k++)
				{
					SNPTABLE *cs = recSNPpos.loc + k;
					if ((cs->chr == chr) && (cs->pos > curpos) && (cs->pos < curpos+(w.ws_col-kScreenPanelMargin*gScreenCnt)/gScreenCnt))
					{
						for (i = 0; i < gScreenCnt; i++)
							printf("\033[%u;%df\033[1;%um\033[1;%um%c\033[0m\n",
								voff,
								(cs->pos-curpos)+i*((w.ws_col-kScreenPanelMargin*gScreenCnt)/gScreenCnt+kColOffset)+kColOffset,
								cs->bgCode,
								cs->fgCode,
								genome[cs->pos-availFromPos-1]); // +kColOffset because of offset line mnumbers +1  -1 because genome is zero base ?
					}
				}
			}
		}


		// print prompts
		printf("\033[%d;1f\033[J",w.ws_row);
		for (i = 0; i < gScreenCnt; i++)
		{
			if (activeScreen == i)
				printf("\033[%d;%df\033[1;47m\033[1;30m%s\033[0m>%.64s",w.ws_row,i*((w.ws_col-kScreenPanelMargin*gScreenCnt)/gScreenCnt+kColOffset)+1,bam[i].displayname,bam[i].info);
			else
			{
				if (!(showmode & kZoomActivePanel))
					printf("\033[%d;%df%s>%.64s",w.ws_row,i*((w.ws_col-kScreenPanelMargin*gScreenCnt)/gScreenCnt+kColOffset)+1,bam[i].displayname,bam[i].info);
			}
		}
		fflush(stdout);

		// process commands
		if (preExecute != NULL)
		{
			c = *preExecute++;
			if (c == 0)
			{
				read(0,&c,1);
				preExecute = NULL;
			}
		}
		else
			read(0,&c,1);
		switch((char)c)
		{
			/* ----- change View mode ----------------------------------------------------------------------------------------- */

			case 'z':
				if (showmode & kZoomActivePanel)
				{
					showmode &= ~kZoomActivePanel;
				}
				else
				{
					showmode |= kZoomActivePanel;
				}
				break;

			case 'c':
				if ((dir[0] != 0) && (bamcnt > 1) && (gScreenCnt < gNbVirtualScreen))
				{
					int newbam = selectbam(bamfiles,bamcnt,bam,activeScreen,w.ws_row-1,w.ws_col / (maxFileNameWidth + kColOffset),maxFileNameWidth+kColOffset,patient);
					showmode |= kSplitScreen;
					sprintf(bam[gScreenCnt].fn,"%s/%s.bam",dir,&bamfiles[newbam*kMaxBamName]);
					if (access(bam[gScreenCnt].fn,R_OK) != 0)
						sprintf(bam[gScreenCnt].fn,"%s/%s",dir,&bamfiles[newbam*kMaxBamName]);
					strcpy(bam[gScreenCnt].ADNIbamname,&bamfiles[newbam*kMaxBamName]);
					bam[gScreenCnt].curbam = newbam;
					err = openBAM(&bam[gScreenCnt],patient);
					if (err)
						goto bail;
					if (readBAM(&bam[gScreenCnt],tag[gScreenCnt],chr,availFromPos) < 0)
						goto bail;
					gScreenCnt++;
				}
				break;

			case 'x':
				if ((dir[0] != 0) && (bamcnt > 1))
				{
					if (showmode & kSplitScreen)
					{
						closeBAM(&bam[activeScreen]);
						for (i=activeScreen+1;i<gScreenCnt;i++)
						{
							memcpy(&bam[i-1],&bam[i],sizeof(BAMFILE));
							memcpy(tag[i-1],tag[i],gMaxTags*sizeof(TAG));
						}
						bam[gScreenCnt-1].in = NULL;
						bam[gScreenCnt-1].idx = NULL;
						bam[gScreenCnt-1].header = NULL;
						gScreenCnt--;
						if (gScreenCnt == 1)
							showmode &= ~kSplitScreen;
						activeScreen = 0;
					}
					// clear screen.
					printf("\033[H\033[J\n");
				}
				break;

			case 'd':
					if (showmode & kDisplayPatientID)
					{
						showmode &= ~kDisplayPatientID;
						for (i=0;i<gScreenCnt;i++)
							DisplayBAMname(&bam[i]);
					}
					else
					{
						showmode |= kDisplayPatientID;
						for (i=0;i<gScreenCnt;i++)
							DisplayPatientID(&bam[i]);
					}
				break;

			case 'e':
				if ((listCN[0] != -1) && (listMCI[0] != -1) && (listAD[0] != -1))
				{
					if (showmode & kExplorePhenotype)
					{
						int i;
						showmode &= ~kExplorePhenotype;
						for (i = 0; i < gScreenCnt; i++)
						{
							closeBAM(&bam[i]);
							bam[i].curbam = i;
						}
						goto awfulHack;
					}
					else
					{
						int i;
						int cn = 0;
						int mci = 0;
						int ad = 0;
						showmode |= kExplorePhenotype;
						for (i = 0; i < gScreenCnt; i++)
						{
							closeBAM(&bam[i]);
							if (i < (gScreenCnt/3))
								bam[i].curbam = listCN[cn++];
							else if (i < 2*(gScreenCnt/3))
								bam[i].curbam = listMCI[mci++];
							else
								bam[i].curbam = listAD[ad++];
						}
						goto awfulHack;
					}
				}
				break;

			case 'f':
				if (showmode & kMaskLowQuality)
					showmode &= ~kMaskLowQuality;
				else
					showmode |= kMaskLowQuality;
				break;

			case '#':
					//tcsetattr(0,TCSANOW,&old);
					if (preExecute != NULL && *preExecute != 0)
						gLowQualityFilterThreshold = *preExecute++;
					else
					{
						printf("\033[%d;1f\033[J\033[1;47m\033[1;30m%s\033[0m",w.ws_row,"Mask mismatches when Quality is <= [type a char]:");
						fflush(stdout);
						gLowQualityFilterThreshold = (char)getchar();
					}
					//tcsetattr(0,TCSANOW,&new);
				break;

			case 'V':
				{
					char input[kMaxFileName];

					bam[0].info[0] = 0;
					if (preExecute != NULL && *preExecute != 0)
					{
						char *s = input;
						int C = *preExecute;
						while (isdigit(C))
						{
							*s++ = C;
							preExecute += 1;
							C = *preExecute;
						}
						if (C == ':' && *(preExecute+1) != 0)
						{
							*s++ = *preExecute++;
							*s++ = *preExecute++;
						}
						*s = 0;
					}
					else
					{
						printf("\033[%d;1f\033[J\033[1;47m\033[1;30m%s\033[0m",w.ws_row,"[Enter position:nt]:");
						fflush(stdout);
						tcsetattr(0,TCSANOW,&old);
						fgets(input,kMaxFileName-1,stdin);
						tcsetattr(0,TCSANOW,&new);
					}
					int nparam = sscanf(input,"%d:%c",&gFilterAllelePos,&gFilterAlleleNT);
					if ((nparam == 2) && (gFilterAllelePos > 0)) // pos.
					{
						showmode |= kFilterAllele;
						gFilterAllelePos -= 1; // make it zero based.
						switch (gFilterAlleleNT)
						{
							case 'a': gFilterAlleleNT = 'A'; break;
							case 't': gFilterAlleleNT = 'T'; break;
							case 'c': gFilterAlleleNT = 'C'; break;
							case 'g': gFilterAlleleNT = 'G'; break;
							case 'A':
							case 'T':
							case 'C':
							case 'G':
							case '-':
							case '*':
								break;
							default: showmode &= ~kFilterAllele; ; break;
						}
						break;
					}
					else
					{
						showmode &= ~kFilterAllele;
						printf("\033[Jinvalid input");
					}
				}
				break;

			case 'C':
				if (showmode & kShowCoverage)
					showmode &= ~kShowCoverage;
				else
					showmode |= kShowCoverage;
				break;

			case 'Z':
				if (showmode & kZoomCoverage)
					showmode &= ~kZoomCoverage;
				else
					showmode |= kZoomCoverage;
				break;

			case 'K':
				if (showmode & kFilterClone)
					showmode &= ~kFilterClone;
				else
					showmode |= kFilterClone;
				goto changechr;
				break;

			case 'O':
				if (showmode & kReadCoverage)
					showmode &= ~kReadCoverage;
				else
					showmode |= kReadCoverage;
				goto changechr;
				break;

			case 'o':
				if (showmode & kShowTagOrdinal)
					showmode &= ~kShowTagOrdinal;
				else
					showmode |= kShowTagOrdinal;
				break;

			case 't':
				if (showmode & kShowTagName)
					showmode &= ~kShowTagName;
				else
					showmode |= kShowTagName;
				break;

			case 'T':
				if (showmode & kHideSameMapping)
					showmode &= ~kHideSameMapping;
				else
					showmode |= kHideSameMapping;
				break;

			case 'q':
				if (showmode & kShowQuality)
					showmode &= ~kShowQuality;
				else
					showmode |= kShowQuality;
				break;

			case 'w':
				if (showmode & kShowSequence)
					showmode &= ~kShowSequence;
				else
					showmode |= kShowSequence;
				break;

			case 'v':
				if (showmode & kShowSNP)
					showmode &= ~kShowSNP;
				else
					showmode |= kShowSNP;
				break;

			/* ----- scroll View ----------------------------------------------------------------------------------------- */
			case 'm':
				{
					int a,b;
					char input[kMaxFileName];

					bam[0].info[0] = 0;
					if (preExecute != NULL && *preExecute != 0)
					{
						char *s = input;
						int C = *preExecute;
						while (isdigit(C))
						{
							*s++ = C;
							preExecute += 1;
							C = *preExecute;
						}
						if (C == ':')
						{
							*s++ = *preExecute++;
							C = *preExecute;
							while (isdigit(C))
							{
								*s++ = C;
								preExecute += 1;
								C = *preExecute;
							}
						}
						*s = 0;
					}
					else
					{
						printf("\033[%d;1f\033[J\033[1;47m\033[1;30m%s\033[0m",w.ws_row,"[Enter new position]:");
						fflush(stdout);
						tcsetattr(0,TCSANOW,&old);
						fgets(input,kMaxFileName-1,stdin);
						tcsetattr(0,TCSANOW,&new);
					}
					int nparam = sscanf(input,"%d:%d",&a,&b);
					if ((nparam == 2) && (a > 0) && (a <= 24)) // chr and pos.
					{
						chr = a;
						curpos = b;
						availFromPos = (b > kHalfGenomicChunk) ? b-kHalfGenomicChunk : 0;
						updateGenomeSeq = 1;
						goto changechr;
					}
					else if (nparam == 1)
					{
						if (scrollView(&availFromPos,&curpos,a-curpos,w.ws_row)) {
							updateGenomeSeq = 1;
							goto changechr;
						}
					}
					else if (gShowAnnotation)
					{
						size_t len = strlen(input);
						char cmd[1024];
						if (len > 0)
							input[len - 1] = 0; // squash newline
						sprintf(cmd,"echo %s/N_NC_*.txt | xargs -n 1 -P 0 grep -H -b -o %s | sed 's+^%s/N_NC_0*++;s+\\.txt++'",gShowAnnotation,input,gShowAnnotation);
						if (verbose)
							fprintf(stderr,"fetching gene name '%s'\n",cmd);
						FILE *f = popen (cmd, "r");
						geneSearchPos.nb = 0;
						geneSearchPos.cur = 0;
						while (fgets(cmd,1024,f) != NULL)
						{
							if (sscanf(cmd,"%d:%d:",&a,&b) == 2)
							{
								if (geneSearchPos.nb == geneSearchPos.max)
								{
									geneSearchPos.max += 1024;
									geneSearchPos.loc = realloc(geneSearchPos.loc, geneSearchPos.max * sizeof(SNPTABLE));
									if (geneSearchPos.loc == NULL)
									{
										fprintf(stderr, "Error: realloc failed\n");
										geneSearchPos.nb = 0;
										geneSearchPos.cur = 0;
										geneSearchPos.max = 0;
										break;
									}
								}
								SNPTABLE *cs = geneSearchPos.loc + geneSearchPos.nb;
								geneSearchPos.nb += 1;
								cs->chr = a;
								cs->pos = b;
								cs->bgCode = 0;
								cs->fgCode = 0;
							}
						}
						pclose(f);
						if (geneSearchPos.nb > 0)
						{
							chr = geneSearchPos.loc[0].chr;
							curpos = geneSearchPos.loc[0].pos;
							availFromPos = (curpos > kHalfGenomicChunk) ? curpos-kHalfGenomicChunk : 0;
							updateGenomeSeq = 1;
							curNPpos = &geneSearchPos;
							goto changechr;
						}
					}
				}
				break;

			/* ----- search genome --------------------------------------------------------------------------------------- */
			case '\\':
				{
					searchGenome[0] = 0;
					genomeSearchPos.nb = 0;
					genomeSearchPos.cur = 0;
					if (preExecute != NULL && *preExecute != 0)
					{
						char *s = searchGenome;
						while (*preExecute != 0 && *preExecute != '\\')
							*s++ = *preExecute++;
						if (*preExecute != 0)
							preExecute += 1;
						*s++ = '\n';
						*s = 0;
					}
					else
					{
						printf("\033[%d;1f\033[J\033[1;47m\033[1;30m%s\033[0m",w.ws_row,"[Enter search regex]:");
						fflush(stdout);
						tcsetattr(0,TCSANOW,&old);
						fgets(searchGenome,512,stdin);
						tcsetattr(0,TCSANOW,&new);
					}
					size_t len = strlen(searchGenome);
					if (len > 0)
					{
						len -= 1;
						searchGenome[len] = 0; // kill newline
					}
					if (len > 0)
					{
						char *s = genome;
						regex_t re;
						regmatch_t rm;
						regcomp(&re, searchGenome, REG_EXTENDED);
						while (regexec(&re, s, 1, &rm, REG_NOTBOL|REG_NOTEOL) == 0)
						{
							if (genomeSearchPos.nb == genomeSearchPos.max)
							{
								genomeSearchPos.max += 1024;
								genomeSearchPos.loc = realloc(genomeSearchPos.loc, genomeSearchPos.max * sizeof(SNPTABLE));
								if (genomeSearchPos.loc == NULL)
								{
									fprintf(stderr, "Error: realloc failed\n");
									genomeSearchPos.nb = 0;
									genomeSearchPos.cur = 0;
									genomeSearchPos.max = 0;
									break;
								}
							}
							SNPTABLE *cs = genomeSearchPos.loc + genomeSearchPos.nb;
							genomeSearchPos.nb += 1;
							cs->chr = chr;
							cs->pos = s - genome + availFromPos + rm.rm_so;
							cs->bgCode = 0;
							cs->fgCode = 0;
							s += rm.rm_eo;
						}
						regfree(&re);
					}
					if (genomeSearchPos.nb > 0) // FIXME - if 0 could use fetchGWI to find a match somewhere
						curNPpos = &genomeSearchPos;
					else
						curNPpos = &SNPpos;
					if (curNPpos->nb > 0)
					{
						int best = find_closest_match(curNPpos,chr,curpos);
						if (best != -1)
							curNPpos->cur = best;
						sprintf(bam[0].info,"%u:%u",curNPpos->cur+1,curNPpos->loc[curNPpos->cur].pos);
						if (curNPpos->loc[curNPpos->cur].chr != chr)
						{
							chr = curNPpos->loc[curNPpos->cur].chr;
							curpos = curNPpos->loc[curNPpos->cur].pos-midscreen;
							availFromPos = (curNPpos->loc[curNPpos->cur].pos > kHalfGenomicChunk) ? curNPpos->loc[curNPpos->cur].pos-kHalfGenomicChunk : 0;
							updateGenomeSeq = 1;
							goto changechr;
						}
						else
						if (scrollView(&availFromPos,&curpos,curNPpos->loc[curNPpos->cur].pos-curpos-midscreen,w.ws_row)) {
							updateGenomeSeq = 1;
							goto changechr;
						}
					}
				}
				break;

			case 'l':
				bam[0].info[0] = 0;
				if (scrollView(&availFromPos,&curpos,midscreen/4,w.ws_row)) {
					updateGenomeSeq = 1;
					goto changechr;
				}
				break;

			case 'L':
				bam[0].info[0] = 0;
				if (scrollView(&availFromPos,&curpos,midscreen,w.ws_row)) {
					updateGenomeSeq = 1;
					goto changechr;
				}
				break;

			case 'j':
				bam[0].info[0] = 0;
				if (scrollView(&availFromPos,&curpos,-midscreen/4,w.ws_row)) {
					updateGenomeSeq = 1;
					goto changechr;
				}
				break;

			case 'J':
				bam[0].info[0] = 0;
				if (scrollView(&availFromPos,&curpos,-midscreen,w.ws_row)) {
					updateGenomeSeq = 1;
					goto changechr;
				}
				break;

			case 'i':
				if (line > 10)
					line-=10;
				else
					line = 1;
				break;

			case 'I':
				if (line > 1)
					line--;
				break;

			case 'k':
				line+=10;
				break;

			/* ----- change SNP ----------------------------------------------------------------------------------------- */
			case 'N':
				if (recSNPpos.nb > 0)  // recorded snps.
				{
						if (recSNPpos.cur < (recSNPpos.nb-1)) recSNPpos.cur++; else recSNPpos.cur = 0;
						SNPTABLE *cs = recSNPpos.loc + recSNPpos.cur;
						sprintf(bam[0].info,"%u:%u",recSNPpos.cur+1,cs->pos);
						int best = find_closest_match(&SNPpos,cs->chr,cs->pos);
						if (best != -1)
							SNPpos.cur = best;
						if (cs->chr != chr)
						{
							chr = cs->chr;
							curpos = cs->pos-midscreen;
							availFromPos = (cs->pos > kHalfGenomicChunk) ? cs->pos-kHalfGenomicChunk : 0;
							updateGenomeSeq = 1;
							line = 1;
							goto changechr;
						}
						else
						if (scrollView(&availFromPos,&curpos,cs->pos-curpos-midscreen,w.ws_row)) {
							updateGenomeSeq = 1;
							line = 1;
							goto changechr;
						}
				}
				break;
			case 'P':
				if (recSNPpos.nb > 0)  // recorded snps.
				{
						if (recSNPpos.cur > 0) recSNPpos.cur--; else recSNPpos.cur = (recSNPpos.nb-1);
						SNPTABLE *cs = recSNPpos.loc + recSNPpos.cur;
						sprintf(bam[0].info,"%u:%u",recSNPpos.cur+1,cs->pos);
						int best = find_closest_match(&SNPpos,cs->chr,cs->pos);
						if (best != -1)
							SNPpos.cur = best;
						if (cs->chr != chr)
						{
							chr = cs->chr;
							curpos = cs->pos-midscreen;
							availFromPos = (cs->pos > kHalfGenomicChunk) ? cs->pos-kHalfGenomicChunk : 0;
							updateGenomeSeq = 1;
							line = 1;
							goto changechr;
						}
						else
						if (scrollView(&availFromPos,&curpos,cs->pos-curpos-midscreen,w.ws_row)) {
							updateGenomeSeq = 1;
							line = 1;
							goto changechr;
						}
				}
				break;

			case 'n':
				if (curNPpos != NULL && curNPpos->nb > 0)  // we have a navigation context
				{
						if (curNPpos->cur < (curNPpos->nb-1)) curNPpos->cur++; else curNPpos->cur = 0;
						sprintf(bam[0].info,"%u:%u",curNPpos->cur+1,curNPpos->loc[curNPpos->cur].pos);
						if (curNPpos->loc[curNPpos->cur].chr != chr)
						{
							chr = curNPpos->loc[curNPpos->cur].chr;
							curpos = curNPpos->loc[curNPpos->cur].pos-midscreen;
							availFromPos = (curNPpos->loc[curNPpos->cur].pos > kHalfGenomicChunk) ? curNPpos->loc[curNPpos->cur].pos-kHalfGenomicChunk : 0;
							updateGenomeSeq = 1;
							line = 1;
							goto changechr;
						}
						else
						if (scrollView(&availFromPos,&curpos,curNPpos->loc[curNPpos->cur].pos-curpos-midscreen,w.ws_row)) {
							updateGenomeSeq = 1;
							line = 1;
							goto changechr;
						}
				}
				break;

			case 'p':
				if (curNPpos != NULL && curNPpos->nb > 0)  // we have a navigation context
				{
						if (curNPpos->cur > 0) curNPpos->cur--; else curNPpos->cur = (curNPpos->nb-1);
						sprintf(bam[0].info,"%u:%u",curNPpos->cur+1,curNPpos->loc[curNPpos->cur].pos);
						if (curNPpos->loc[curNPpos->cur].chr != chr)
						{
							chr = curNPpos->loc[curNPpos->cur].chr;
							curpos = curNPpos->loc[curNPpos->cur].pos-midscreen;
							availFromPos = (curNPpos->loc[curNPpos->cur].pos > kHalfGenomicChunk) ? curNPpos->loc[curNPpos->cur].pos-kHalfGenomicChunk : 0;
							updateGenomeSeq = 1;
							line = 1;
							goto changechr;
						}
						else
						if (scrollView(&availFromPos,&curpos,curNPpos->loc[curNPpos->cur].pos-curpos-midscreen,w.ws_row)) {
							updateGenomeSeq = 1;
							line = 1;
							goto changechr;
						}
				}
				break;

			case 'r':
				if (curNPpos != NULL && curNPpos->nb > 0)  // we have a navigation context
				{
					if (recSNPpos.nb == recSNPpos.max)
					{
						recSNPpos.max += 1024;
						recSNPpos.loc = realloc(recSNPpos.loc, recSNPpos.max * sizeof(SNPTABLE));
						if (recSNPpos.loc == NULL)
						{
							fprintf(stderr, "Error: realloc failed\n");
							break;
						}
					}
					SNPTABLE *cs = recSNPpos.loc + recSNPpos.nb;
					recSNPpos.nb += 1;
					*cs = curNPpos->loc[curNPpos->cur];
					cs->bgCode = 45;
					cs->fgCode = 37; // white
				}
				break;

			case 'R':
				{
					unsigned int i;
					fprintf(stderr,"chr\tchrPos\n");
					for (i = 0; i < recSNPpos.nb; i++)
						fprintf(stderr,"chr%d\t%d\t%u\t%u\n",recSNPpos.loc[i].chr,recSNPpos.loc[i].pos,recSNPpos.loc[i].bgCode,recSNPpos.loc[i].fgCode);
				}
				break;

			case 'g':
				{
					int a,b;
					char input[kMaxFileName];

					bam[0].info[0] = 0;
					if (preExecute != NULL && *preExecute != 0)
					{
						char *s = input;
						int C = *preExecute;
						if (C == 'G' || C == 'M' || C == 'S')
						{
							*s++ = C;
							preExecute += 1;
							C = *preExecute;
						}
						while (isdigit(C))
						{
							*s++ = C;
							preExecute += 1;
							C = *preExecute;
						}
						if (C == ':')
						{
							*s++ = *preExecute++;
							C = *preExecute;
							while (isdigit(C))
							{
								*s++ = C;
								preExecute += 1;
								C = *preExecute;
							}
						}
						*s = 0;
					}
					else
					{
						printf("\033[%d;1f\033[J\033[1;47m\033[1;30m%s%d%s\033[0m",w.ws_row,"[Enter item number (1..",curNPpos->nb,") or item position]:");
						fflush(stdout);
						tcsetattr(0,TCSANOW,&old);
						fgets(input,kMaxFileName-1,stdin);
						tcsetattr(0,TCSANOW,&new);
					}
					int nparam;
					if (input[0] == 'G' || input[0] == 'M' || input[0] == 'S')
					{
						nparam = sscanf(input + 1,"%d:%d",&a,&b);
						switch (input[0])
						{
							case 'G' : curNPpos = &geneSearchPos; break;
							case 'M' : curNPpos = &genomeSearchPos; break;
							case 'S' : curNPpos = &SNPpos; break;
						}
					}
					else
						nparam = sscanf(input,"%d:%d",&a,&b);
					if (nparam == 2) // chr and pos.
					{
						int best = find_closest_match(curNPpos,a,b);
						if (best != -1)
							curNPpos->cur = best;
					}
					else if ((nparam == 1) && (a > 0) && (a <= curNPpos->nb))
					{
						curNPpos->cur = a - 1;
					}
					else
					{
						printf("\033[JValid index are between 1 and %d", curNPpos->nb);
					}
					sprintf(bam[0].info,"%u:%u",curNPpos->cur+1,curNPpos->loc[curNPpos->cur].pos);
					if (curNPpos->loc[curNPpos->cur].chr != chr)
					{
						chr = curNPpos->loc[curNPpos->cur].chr;
						curpos = curNPpos->loc[curNPpos->cur].pos-midscreen;
						availFromPos = (curNPpos->loc[curNPpos->cur].pos > kHalfGenomicChunk) ? curNPpos->loc[curNPpos->cur].pos-kHalfGenomicChunk : 0;
						updateGenomeSeq = 1;
						goto changechr;
					}
					else
					if (scrollView(&availFromPos,&curpos,curNPpos->loc[curNPpos->cur].pos-curpos-midscreen,w.ws_row)) {
						updateGenomeSeq = 1;
						goto changechr;
					}
				}
				break;
			/* ----- change bam ----------------------------------------------------------------------------------------- */
			case 'a':
				if (bamcnt > 1)
				{
					int newbam;
					if (bam[activeScreen].curbam > 0)
						newbam = bam[activeScreen].curbam-1;
					else
						newbam = bamcnt - 1;
					if (changeActiveScreenBam(&bam[activeScreen],tag[activeScreen],newbam,bamfiles,dir,chr,availFromPos,patient) != 0)
						goto bail;
				}
				break;

			case 's':
				if (bamcnt > 1)
				{
					int newbam;
					if (bam[activeScreen].curbam < (bamcnt-1))
						newbam = bam[activeScreen].curbam+1;
					else
						newbam = 0;
					if (changeActiveScreenBam(&bam[activeScreen],tag[activeScreen],newbam,bamfiles,dir,chr,availFromPos,patient) != 0)
						goto bail;
				}
				break;

			case 'A':
				if (bamcnt > 1)
				{
					EXECUTIONPLAN ep[gNbVirtualScreen];
					pthread_t thread[gNbVirtualScreen];
					int lastbam = bamcnt-1;
					for (i = 0; i < gScreenCnt; i++)
					{
						closeBAM(&bam[i]);
						if (showmode & kExplorePhenotype)
						{

							if (i < (gScreenCnt/3))
								bam[i].curbam = getPreviousPatientKind(listCN,bam[i].curbam,(gScreenCnt/3));
							else if (i < 2*(gScreenCnt/3))
								bam[i].curbam = getPreviousPatientKind(listMCI,bam[i].curbam,(gScreenCnt/3));
							else
								bam[i].curbam = getPreviousPatientKind(listAD,bam[i].curbam,(gScreenCnt/3));
						}
						else
						{
							bam[i].curbam -= gScreenCnt;
							if (bam[i].curbam < 0)
								bam[i].curbam = lastbam--;
						}
						ep[i].bam = &bam[i];
						ep[i].tag = tag[i];
						ep[i].chr = chr;
						ep[i].pos = availFromPos;
						sprintf(bam[i].fn,"%s/%s.bam",dir,&bamfiles[bam[i].curbam*kMaxBamName]);
						if (access(bam[i].fn,R_OK) != 0)
							sprintf(bam[i].fn,"%s/%s",dir,&bamfiles[bam[i].curbam*kMaxBamName]);
						strcpy(bam[i].ADNIbamname,&bamfiles[bam[i].curbam*kMaxBamName]);
						err = openBAM(&bam[i],patient);
						if (err)
							goto bail;
						if (pthread_create (&thread[i], NULL, &readBAMthread, &ep[i]) == 0) {
							if (verbose > 1) { fprintf(stderr,"Launching thread %u with %s\n",i,bam[i].fn) ;}
						}
					}

					for (i = 0; i < gScreenCnt; i++)
						if (pthread_join (thread[i], NULL) == 0) { if (verbose > 1) { fprintf(stderr,"Finished thread %u\n",i) ;} }

					for (i = 0; i < gScreenCnt; i++)
					{
						if (ep[i].rslt < 0)
							goto bail;
					}
					compareAlignments(tag[0],tag[1],bam[0].retained,bam[1].retained);

				}
				break;

			case 'S':
				if (bamcnt > 1)
				{
					EXECUTIONPLAN ep[kMaxVirutalScreen];
					pthread_t thread[kMaxVirutalScreen];
					int firstbam = 0;
					for (i = 0; i < gScreenCnt; i++)
					{
						closeBAM(&bam[i]);
						if (showmode & kExplorePhenotype)
						{

							if (i < (gScreenCnt/3))
								bam[i].curbam = getNextPatientKind(listCN,bam[i].curbam,(gScreenCnt/3));
							else if (i < 2*(gScreenCnt/3))
								bam[i].curbam = getNextPatientKind(listMCI,bam[i].curbam,(gScreenCnt/3));
							else
								bam[i].curbam = getNextPatientKind(listAD,bam[i].curbam,(gScreenCnt/3));
						}
						else
						{
							bam[i].curbam+=gScreenCnt;
							if (bam[i].curbam >= bamcnt)
								bam[i].curbam = firstbam++;
						}
					}
		awfulHack:;

					for (i = 0; i < gScreenCnt; i++)
					{
						ep[i].bam = &bam[i];
						ep[i].tag = tag[i];
						ep[i].chr = chr;
						ep[i].pos = availFromPos;
						sprintf(bam[i].fn,"%s/%s.bam",dir,&bamfiles[bam[i].curbam*kMaxBamName]);
						if (access(bam[i].fn,R_OK) != 0)
							sprintf(bam[i].fn,"%s/%s",dir,&bamfiles[bam[i].curbam*kMaxBamName]);
						strcpy(bam[i].ADNIbamname,&bamfiles[bam[i].curbam*kMaxBamName]);
						err = openBAM(&bam[i],patient);
						if (err)
							goto bail;
						if (pthread_create (&thread[i], NULL, &readBAMthread, &ep[i]) == 0) {
							if (verbose > 1) { fprintf(stderr,"Launching thread %u with %s\n",i,bam[i].fn) ;}
						}
					}

//					for (i = 0; i < 2; i++)
					for (i = 0; i < gScreenCnt; i++)
						if (pthread_join (thread[i], NULL) == 0) { if (verbose > 1) { fprintf(stderr,"Finished thread %u\n",i) ;} }

//					for (i = 0; i < 2; i++)
					for (i = 0; i < gScreenCnt; i++)
					{
						if (ep[i].rslt < 0)
							goto bail;
					}
					compareAlignments(tag[0],tag[1],bam[0].retained,bam[1].retained);

				}
				break;

			case '/':
				if (bamcnt > 1)
				{
					int newbam = selectbam(bamfiles,bamcnt,bam,activeScreen,w.ws_row-1,w.ws_col / (maxFileNameWidth + kColOffset),maxFileNameWidth+kColOffset,patient);
					if (newbam != bam[activeScreen].curbam)
					{
						if (changeActiveScreenBam(&bam[activeScreen],tag[activeScreen],newbam,bamfiles,dir,chr,availFromPos,patient) != 0)
							goto bail;
					}
				}
				break;

			case 'U':
				if (alignerConfigFile[0])
				{
					err = readGTL(&bam[0],tag[0],chr,curpos,1);
					if (err < 0)
					{
						fprintf(stderr,"EXIT:chr%d\n",err);
						goto bail;
					}
				}
				break;


			case 'E':
			case 'F':
				{
					int line;
					char input[kMaxFileName];
					if (preExecute != NULL && *preExecute != 0)
					{
						char *s = input;
						int C = *preExecute;
						while (isdigit(C))
						{
							*s++ = C;
							preExecute += 1;
							C = *preExecute;
						}
						*s = 0;
					}
					else
					{
						printf("\033[%d;1f\033[J\033[1;47m\033[1;30m%s\033[0m",w.ws_row,"[Enter Line number:");
						fflush(stdout);
						tcsetattr(0,TCSANOW,&old);
						fgets(input,kMaxFileName-1,stdin);
						tcsetattr(0,TCSANOW,&new);
					}
					sscanf(input,"%d",&line);
					// FIXME - this fails when there are indels
					if ((char)c =='E')
						printtagfromline(tag[activeScreen],bam[activeScreen].retained,chr,line);
					else
					{
						if (line != 0)
							fileappendtagfromline(tag[activeScreen],bam[activeScreen].retained,chr,line,1);
						else
						{
							int i;
							for (i=1; i<bam[activeScreen].retained; i++)
								fileappendtagfromline(tag[activeScreen],bam[activeScreen].retained,chr,i,0);
						}
					}
				}
			break;

			/* ----- change active screen panel ---------------------------------------------------------------------------- */

			case (char)0x09:
				activeScreen ++;
				if (activeScreen >= gScreenCnt)
					activeScreen = 0;
				break;

			case '&':
				if (curNPpos != NULL && curNPpos->nb > 0)  // we have a navigation context
				{
					int as;
					int maxlcnt = 0;
					for (as=0; as<gScreenCnt; as++)
					{
						printchunk(tag[as],bam + as,chr,line,availFromPos,curpos,(w.ws_col-kScreenPanelMargin*1)/gScreenCnt,  kVirtScreenMaxLines-1  ,&screenpanel[as],bam[as].color);
						if (screenpanel[as].printed > maxlcnt)
							maxlcnt = screenpanel[as].printed;
					}
					for (as=0; as<gScreenCnt; as++)
					{
						char vsfn[kMaxFileName];
						int idx;

						sprintf(vsfn,"/tmp/c%dp%d_%s.txt",curNPpos->loc[curNPpos->cur].chr,curNPpos->loc[curNPpos->cur].pos,bam[as].displayname);
						flushscreenPanelToFile(vsfn,screenpanel,as,maxlcnt);

						sprintf(vsfn,"/tmp/c%dp%d_%s.tsv",curNPpos->loc[curNPpos->cur].chr,curNPpos->loc[curNPpos->cur].pos,bam[as].displayname);
						FILE *of = fopen(vsfn,"w");
						fprintf(of,"chr\tpos\tnt\tA\tC\tG\tT\tdel\tins\n");
						for(idx = 0; idx < (w.ws_col-kScreenPanelMargin*1)/gScreenCnt; idx++)
						{
							unsigned int p = curpos - bam[as].startPos + idx;
							unsigned int cp = kCoverageElements * p;
							fprintf(of,"%d\t%d\t%c\t%hu\t%hu\t%hu\t%hu\t%hu\t%hu\n",curNPpos->loc[curNPpos->cur].chr,curpos+idx+1,genome[curpos-availFromPos+idx],bam[as].coverage[cp+0],bam[as].coverage[cp+1],bam[as].coverage[cp+2],bam[as].coverage[cp+3],bam[as].coverage[cp+4],bam[as].coverage[cp+5]);
						}
						fclose(of);
					}
				}
				break;
			
			/* -----  Help -------------------------------------------------------------------------------------------------- */
			case '?':
			case 'h':
				help();
				break;

			/* -----  Quit -------------------------------------------------------------------------------------------------- */
			case 'Q':
				// clear screen.
				//printf("\033[H\033[J");
				// or not
				printf("\n");
				goto bail;


autoprocess:;
			case '@':
				if (linkedSNPfn[0] != 0)	{

					int jjjj = 0;
					FILE *f = fopen(linkedSNPfn,"w");
					//cut -c10,83,88- /tmp/vscreen_dump.txt  | grep -v "\."  | less
					while(1)
					{
						VIRTUALSCREEN vs;
						err += posix_memalign((void **)&vs.virtualscreen, 64, (12*512*kVirtScreenMaxLines)* sizeof(char)); // some space for color chars.
						if(err == 0)
						{
							if (f)
							{
								//scrollView(&availFromPos,&curpos,45415635-curpos,w.ws_row);
								for (i = 0; i < gScreenCnt; i++)
								{
									printchunkBW(tag[i],bam[i].retained,chr,0,availFromPos,linkedSNP1-1,(linkedSNP2-linkedSNP1+3),(kVirtScreenMaxLines-1),&vs);
									dumpVirtualScreenToFile(&vs,f,bam[i].displayname);
								}
								free(vs.virtualscreen);
							}
						}

						if (jjjj++ == 20)
						//if (jjjj++ == 5)
							break;
						// this section is duplicated from case 'S'
						printf("Reading BAM files #%d\n",jjjj);
						if (bamcnt > 1)
						{
							EXECUTIONPLAN ep[gNbVirtualScreen];
							pthread_t thread[gNbVirtualScreen];
							int firstbam = 0;
							for (i = 0; i < gScreenCnt; i++)
							{
								closeBAM(&bam[i]);
								if (showmode & kExplorePhenotype)
								{

									if (i < (gScreenCnt/3))
										bam[i].curbam = getNextPatientKind(listCN,bam[i].curbam,(gScreenCnt/3));
									else if (i < 2*(gScreenCnt/3))
										bam[i].curbam = getNextPatientKind(listMCI,bam[i].curbam,(gScreenCnt/3));
									else
										bam[i].curbam = getNextPatientKind(listAD,bam[i].curbam,(gScreenCnt/3));
								}
								else
								{
									bam[i].curbam+=gScreenCnt;
									if (bam[i].curbam >= bamcnt)
										bam[i].curbam = firstbam++;
								}
							}

							for (i = 0; i < gScreenCnt; i++)
							{
								ep[i].bam = &bam[i];
								ep[i].tag = tag[i];
								ep[i].chr = chr;
								ep[i].pos = availFromPos;
								sprintf(bam[i].fn,"%s/%s.bam",dir,&bamfiles[bam[i].curbam*kMaxBamName]);
								if (access(bam[i].fn,R_OK) != 0)
									sprintf(bam[i].fn,"%s/%s",dir,&bamfiles[bam[i].curbam*kMaxBamName]);
								strcpy(bam[i].ADNIbamname,&bamfiles[bam[i].curbam*kMaxBamName]);
								err = openBAM(&bam[i],patient);
								if (err)
									goto bail;
								if (pthread_create (&thread[i], NULL, &readBAMthread, &ep[i]) == 0) {
									if (verbose > 1) { fprintf(stderr,"Launching thread %u with %s\n",i,bam[i].fn) ;}
								}
							}

							for (i = 0; i < gScreenCnt; i++)
								if (pthread_join (thread[i], NULL) == 0) { if (verbose > 1) { fprintf(stderr,"Finished thread %u\n",i) ;} }

							for (i = 0; i < gScreenCnt; i++)
							{
								if (ep[i].rslt < 0)
									goto bail;
							}
							compareAlignments(tag[0],tag[1],bam[0].retained,bam[1].retained);

						}

					} // jjjj
					fclose(f);
					goto bail;


				}
				break;

		} // switch

	} // while(keeprunning)


// --------------------------------------------------------------------------------------------------
//                All things have one end.
// --------------------------------------------------------------------------------------------------

bail:
	tcsetattr(0,TCSANOW,&old);

	// close files, free and return
	for (i = 0; i < gScreenCnt; i++)
	{
		closeBAM(&bam[i]);
		if (tag[i])
			free(tag[i]);
	}
	for (i = 0; i < gNbVirtualScreen; i++)
	{
		free(screenpanel[i].virtualscreen);
	}

	if (bamfiles)
		free(bamfiles);

	free(SNPpos.loc);
	free(recSNPpos.loc);
	free(geneSearchPos.loc);
	free(genomeSearchPos.loc);

#ifdef DEBUGOUTPUT
	fclose(debug);
#endif

	return err;
}

// --------------------------------------------------------------------------------------------------
static void usage(void)
{
#ifdef ROBIN
		printf("usage:\n\n");
		printf("samtools ADVIEW -D width [-S] -1 InputFile  -c chromosome -p position [-q] [-m maxtag]\n\n");
		printf("           -D width              : dump width of full alignment of current position to stdout instead of using interactive viewer\n");
		printf("           -S                     : BAM file contains single-end reads\n");
		printf("           -1 InputFile           : BAM file (sorted, and indexed)\n");
		printf("           -c chromosome          : initial chromosome to show [1..24]\n");
		printf("           -p position            : postion on chromosome\n");
		printf("           -q                     : write quality instead of nucleotides\n");
		printf("           -m maxtag              : increase for big files (default = 524288)\n");
		printf("           -l num                 : print every num line\n");
		printf("\nREMARK: bam file must be sorted and indexed.\n");
		printf("to index, use the following command   samtools index file.bam");
		printf("\n");
		printf("\nexamples:\n");
		printf("./samtools ADVIEW -1 IonXpress_017_rawlib.bam -c 5 -p 112173810  -S -D\n");
		printf("Authors: Nicolas Guex and Christian Iseli, 2014-2019\nThis program comes with ABSOLUTELY NO WARRANTY.\nThis is free software, released under GPL2+ and you are welcome to redistribute it under certain conditions.\n");
		printf("CONTACT: Nicolas.Guex@unil.ch  and  Christian.Iseli@unil.ch \n");
		printf("Version 1.007\n\n");
		printf("note: making use of samtools\n");
		printf("Copyright (C) 2008-2019 Genome Research Ltd.\n");
		printf("Portions copyright (C) 2009, 2011, 2012 Broad Institute.\n");
#else
		printf("usage:\n\n");
		printf("samtools ADVIEW  [-h hightlight_file] ([-d directory [-n screens] | [-1 InputFile [-2 InputFile [-3 InputFile etc...]]]) -c chromosome -p position [ -s snpfile ]  [-v level] [-A position -B position -L outputfile]\n\n");
		printf("           -C                     : do not try to compare read names in dual screen mode\n");
		printf("                                  : indicates/instruct the viewer to fetch patient id, and SNPs using Alex rest API.\n");
		printf("           -h highlight_file      : load a tab delimited text file with two columns patient_id and color\n");
		printf("                                    this controles the highlight of the display header 0=no_highlight (default)\n");
		printf("                                    1=green; 2=yellow; 3=red\n");
		printf("           -d Directory           : containing BAM files (sorted, and indexed)\n");
		printf("           -1..9 InputFile        : BAM file (sorted, and indexed)\n");
		printf("           -c chromosome          : initial chromosome to show [1..24]\n");
		printf("           -p position            : postion on chromosome\n");
		printf("           -s snpfile             : file containing location of snps to visualize\n");
		printf("                                    when present it supercedes the http connection to fetch Alex SNP data.\n");
		printf("                                    fist line considered the header and is ignored \n");
		printf("                                    should be tab delimited  chr position\n");
		printf("                                    chromosomes should be chr1..chr22,chrX or chrY\n");
		printf("                                    and yes, it could be used as bookmarks ;-)\n");
		printf("           -n screens             : specifies the initial number of split screens [1..%u]\n",kMaxVirutalScreen);
		printf("           -G GTL genome          : specifies the genome file GTLdecompress and GTLfetch will use (default: %s)\n",gGTLgenome);
		printf("           -g genome              : specifies the genome file Chris' fetch will use (default is to use GTLfetch)\n");
		printf("           -v level               : specifies the verbose level; default is 0\n\n");
		printf("           -A position            : specifies something to do with linked SNPs...\n\n");
		printf("           -B position            : specifies something to do with linked SNPs...\n\n");
		printf("           -L file                : specifies linked SNPs file name (?)\n\n");
		printf("           -S                     : BAM file contains single-end reads\n\n");
		printf("           -D                     : dump full alignment of current position to stdout instead of using interactive viewer\n");
		printf("           -z                     : no 'chr' in BAM names\n");
		printf("           -a dir                 : load and show annotations\n");
		printf("           -U file                : config file for real time realignment using Thierry's aligner\n");
		printf("           -e commands            : run all commands specified in the string as if they were typed by the user, before reading user input\n");
		printf("\nexamples:\n");
		printf("./samtools ADVIEW -1 data.bam -c 19 -p 45409011\n");
		printf("\n");
		printf("Authors: Nicolas Guex and Christian Iseli, 2014-2019\nThis program comes with ABSOLUTELY NO WARRANTY.\nThis is free software, released under GPL2+ and you are welcome to redistribute it under certain conditions.\n");
		printf("CONTACT: Nicolas.Guex@unil.ch  and  Christian.Iseli@unil.ch \n");
		printf("Version 1.007\n\n");
		printf("note: making use of samtools\n");
		printf("Copyright (C) 2008-2019 Genome Research Ltd.\n");
		printf("Portions copyright (C) 2009, 2011, 2012 Broad Institute.\n");
#endif
}
// --------------------------------------------------------------------------------------------------
static void help(void)
{
	int c;

		printf("\033[H\033[J");
		fflush(stdout);

		printf("\nCommands:\n");
		printf("h     display this help\n");
		printf("c     create one additional screen, user need to select a BAM (need option -d)\n");
		printf("x     close the active screen (if more than one bam shown; need option -d)\n");
		printf("z     zoom/unzoom active screen panel\n");
		printf("<tab> activate next screen panel (name shown in inverse video)\n");
		printf("Q   quit\n");
		printf("\n");
		printf(" i    scroll view up by 10 rows\n");
		printf("j l   scroll view to the left/right by 10 nt\n");
		printf("J L   scroll view to the left/right by 100 nt\n");
		printf(" k    scroll view down 10 lines\n");
		printf(" m    move view (wait for a new   position or   chromosome:position)\n");
		printf("      example:    19:45409011\n");
		printf("\n");
		printf("C   coverage on 1/3 of screen height; max is 2 times median of the loaded chunk (view / hide)\n");
		printf("K   toggle filtering clones in the coverage view\n");
		printf("Z   zoom on coverage (described in the line above)\n");
		printf("O   toggle coverage on Overall read pair (coverage centric) or per read (variation centric) (see C and Z above)\n");
		printf("q   quality (view / hide)\n");
		printf("v   variant SNP (view / hide)\n");
		printf("V   variant allele (filter view to show only tags with specific NT at a given position)\n");
		printf("o   tag ordinal (view / hide)\n");
		printf("t   tag name (view / hide)\n");
		printf("T   when tagnames displayed, and only two files loaded, will not display tagpairs mapped at same place\n");
		printf("d   display toggle between patient ID and BAM filename\n");
		printf("f   filter nucleotides with low quality (#)\n");
		printf("#   set up the quality score value used to mask nucleotides with low quality\n");
		printf("w   actual sequence (view / hide)\n");
		printf("\\  highlights a sequence fragment in the genome\n");
		printf("U   update alignment of second view panel (soft. needs to be launched with -U file)\n");
//		printf("t   toggle view between normal and compact modes\n");
//		printf("    in compact mode, each printed line can have at most 2 mismatches\n");
//		printf("    (check is done only on the visible portion of each tag)\n");
//		printf("    in compact mode, only paired-end tags for which at least one of the pair\n");
//		printf("    is a perfect match are shown.\n");
		printf("\n");
		printf("n   center view on next     SNP  of supplied -s file\n");
		printf("p   center view on previous SNP  of supplied -s file\n");
		printf("N   center view on next     SNP  of the active screen panel\n");
		printf("P   center view on previous SNP  of the active screen panel\n");
		printf("g   select a SNP (wait for a SNP number or search a SNP by chromosome:position)\n");
		printf("r   record current SNP position\n");
		printf("R   output list of recorded SNPs to stderr\n");
		printf("\n");
		printf("e   turn on/off the Exploration mode. (arrange screenPanels/3 patients of each category for browsing)\n");
		printf("a   switch active screen to the previous patient (requires use of the option -d at launch)\n");
		printf("A   switch all   screens to the previous patient (requires use of the option -d at launch)\n");
		printf("s   switch active screen to the next patient (requires use of the option -d at launch)\n");
		printf("S   switch all   screens to the next patient (requires use of the option -d at launch)\n");
		printf("/   select bam file from directory (requires use of the option -d at launch)\n");
		printf("E   Extract the tag pair from a given line of the active screen as fasta\n");
		printf("F   FastQ append tag pair from a given line to ./dbg_R1.fq and ./dbg_R2.fq (use 0 for all lines)\n");
		printf("&   write content of current display into files, 1 file per sample, in /tmp\n");
		printf("\n  press any key to exit this screen\n");

		read(0,&c,1);

}
// vim: tabstop=2 shiftwidth=2
// --------------------------------------------------------------------------------------------------
