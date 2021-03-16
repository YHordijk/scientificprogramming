/*
 *     Copyright (c) 2019 by the authors of DIRAC.
 *     All Rights Reserved.
 *
 *     This source code is part of the DIRAC program package.
 *     It is provided under a written license and may be used,
 *     copied, transmitted, or stored only in accordance to the
 *     conditions of that written license.
 *
 *     In particular, no part of the source code or compiled modules may
 *     be distributed outside the research group of the license holder.
 *     This means also that persons (e.g. post-docs) leaving the research
 *     group of the license holder may not take any part of Dirac,
 *     including modified files, with him/her, unless that person has
 *     obtained his/her own license.
 *
 *     For information on how to get a license, as well as the
 *     author list and the complete list of contributors to the
 *     DIRAC program, see: http://www.diracprogram.org
 */

/* -*- Mode: C; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/*
 * Implementation of Embedding Functionals
 *
 * Written by
 * Christoph Jacob, jacob@few.vu.nl, 2004
 *
 */


/* FIXME: portability issues with the str... functions. Some libc
   provide them only for BSD compatibility */
#define _DEFAULT_SOURCE 1
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#define __CVERSION__
#include "general.h"
#include "functionals.h"

/* 
 * Embedding functional definition:
 *
 * An "embedding functional" consists of an exchange-correlation functional
 * and a kinetic energy functional and depends on TWO densities:
 * the density of a part of the system (rho1 - the "interesting" part) and 
 * the total density (rho_tot = rho1 + rho2 - the complete system)
 *
 * It REPLACES the "normal" exchange-correlation functional for a 
 * calculation, in which rho1 is varying (to be calculated) and rho_tot 
 * is frozen. 
 *
 * Energy terms are calculated for the WHOLE system.
 * 
 * The DFT-potential (first order) REPLACES the "normal" exchange 
 * correlation potential. It is intended for solving the KS-equations
 * in region I, with a frozen rho_tot (or rho_2, for the potential it
 * does not matter).
 * 
 * I havent thougt about the second and third order terms yet.
 */

/* 
 * "Embedding Interaction" functional definition:
 *
 * An "Embedding Interaction" functional is quite similar to the 
 * embedding functional described above. The only difference is,
 * that it only calculates the INTERACTION contributions, not the
 * total energy/potential terms.
 *
 * The calculated interaction contributions have to be added to
 * the contributions from region I and from region II. For the
 * interaction part it does not matter, whether rho2 or rho_total is 
 * frozen, this has to be taken care of when calculating the
 * contributions from region I and II
 *
 * Energy terms are the INTERACTION ENERGY only. E(region I) and
 * E(region II) must be added to get the total energy. (This functional
 * only takes care of the XC and kinetic energy parts, no electrostatic
 * interactions are included!) 
 * 
 * The DFT-potential (fisrt order) has to be added to the "normal" 
 * exchange correlation potential. 
 * It is intended for solving the KS-equations in region I, with 
 * a frozen rho_tot (or rho_2, for the potential it does not matter).
 * 
 * As for the embedding, I havent thougt about the second and third 
 * order terms yet.
 */
extern void quit_(const char* str);
static Functional   *embedding_xc_functional = &LDAFunctional;          /* the XC-functional  */
static Functional   *embedding_kin_functional = &TF_KinFunctional;      /* the kinetic energy functional */

/*
 * EmbeddingFunctional and emb_* functions:
 *
 * This is a generic embedding functional in a "Functional" structure,
 * that has to be selected in the DFT part when embedding is done.
 *
 * It provides all the things a normal functional can do, only that the
 * density passed in "dp" has to provide "subsystem_dp"
 * 
 */

static int
emb_isgga (void)
{
	return (embedding_xc_functional->is_gga () || embedding_kin_functional->is_gga () );
}


/* 
 * Reads configuration line for embedding potential:
 * Format of this line:
 *   KIN=[kinetic energy function name, optional parameters] XC=[XC-functional name, optional parameters]
 *   (this follows the name "embedding" in **HAMILTONIAN, *DFT, .DFTFUN)
 *
 * FIXME: this is terribly ugly and way to much work (but the whole Dirac-input is ...) 
 */
static int
emb_read (const char *conf_line)
{
	char *line;                                 /* a copy of the conf_line */
	char *kin_line = NULL, *xc_line = NULL;     /* the KIN=... and XC=... lines, without prefix */
	const char *kin_pos, *xc_pos;               /* pointers to "KIN=" and "XC=" */

	char *pos;
	int success = 1;                            /* our return value */
	int found;
	int i, len;

	/* start processing the conf_line */
	line = strdup (conf_line);
	if (!line)
		return 1;

	/* input is case insensitive, so turn conf_line to uppercase */
	for (pos = line; *pos; pos++) 
		*pos = (char) toupper ((int) *pos);

	/* now first find "KIN=" and "XC=" in conf_line */
	kin_pos = strstr(line, "KIN");
	if (kin_pos && !(kin_pos[3] == '=') )
		kin_pos = NULL;

	xc_pos = strstr(line, "XC");
	if (xc_pos && !(xc_pos[2] == '=') )
		xc_pos = NULL;

	/* now see what we've got and 
	   save the KIN= .... and XC = .... strings (kin_line and xc_line) */
	
	if ( kin_pos && xc_pos) {
		if (kin_pos < xc_pos) {
			xc_line = strdup (xc_pos+3);
			len = xc_pos - (kin_pos+4);
			kin_line = malloc (len+1);
			strncpy (kin_line, kin_pos+4, len);
			kin_line[len] = '\0';
		} else {
			kin_line = strdup (kin_pos+4);
			len = kin_pos - (xc_pos+3);
			xc_line = malloc (len+1);
			strncpy (xc_line, xc_pos+3, len);
			xc_line[len] = '\0';
		}
	} else if ( kin_pos && !xc_pos ) {
		kin_line = strdup (kin_pos+4); 
	} else if ( xc_pos && !kin_pos ) {
		xc_line = strdup (xc_pos+3);
	} else {
		kin_line = xc_line = NULL;
	}

	/* and finally process kin_line and xc_line (if they exist) */

	if (kin_line) {
		for (pos = kin_line; isspace ((int) *pos); pos++);

		found = 0;
		for (i = 0; !found && available_kin_functionals[i]; i++) {
			len =  strlen (available_kin_functionals[i]->name);
			if ( !strncasecmp (available_kin_functionals[i]->name, pos, len) ) {    /* whole name matches */
				found = 1;
				embedding_kin_functional = available_kin_functionals[i];
				success = embedding_kin_functional->read ? 
					                   embedding_kin_functional->read (kin_line + len) : 1;
			} else if ( !strncasecmp (available_kin_functionals[i]->name+4, pos, len-4) ) {  
				found = 1;
				embedding_kin_functional = available_kin_functionals[i];
				success = embedding_kin_functional->read ? 
					                   embedding_kin_functional->read (kin_line +4 +len) : 1;
			}

		}
		if (!found) {
			fort_print ("Unknown kinetic energy funtional '%s'. Aborting. \n", pos);
			success = 0;
		}

		free (kin_line);
	}

	if (success && xc_line) {
		for (pos = xc_line; isspace ((int) *pos); pos++);

		found = 0;
		for (i = 0; !found && available_functionals[i]; i++) {
			if ( !strncasecmp (available_functionals[i]->name, pos, 
				       strlen (available_functionals[i]->name) ) ) {
				found = 1;
				embedding_xc_functional = available_functionals[i];
				success = embedding_xc_functional->read ?
					embedding_xc_functional->read (xc_line + 
									strlen(embedding_xc_functional->name)) : 1;
			}
		}
		if (!found) {
			fort_print ("Unknown exchange-correlation funtional '%s'. Aborting. \n", pos);
			success = 0;
		}

		free (xc_line);
	}

	free (line);
	return success;
}

static void
emb_report (void)
{
	fort_print ("\n     Embedding Functional:");
	fort_print ("          Exchange-correlation functional: %20s", embedding_xc_functional->name);
	if (embedding_xc_functional->report) 
		embedding_xc_functional->report ();
	fort_print ("          Kinetic energy functional:       %20s \n", embedding_kin_functional->name);
	if (embedding_kin_functional->report)
		embedding_kin_functional->report ();
}

static real
emb_energy (const FunDensProp *dp)
{
	/* the total density is frozen, so we know it and the corresponding KS-orbitals 
	   The total energy is just the DFT-energy of rho_total, we dont need any embedding here! */

	return embedding_xc_functional->func (dp); 			
}

static void
emb_first (FunFirstFuncDrv *ds, real factor, const FunDensProp *dp)
{
	if (!dp->subsystem[0]) {
		fort_print ("emb_first called but dp->subsystem[0] not initialized\n");
		fort_print ("This should not happen ... \n");
		quit_ ("Error in emb_first");
	} else {
		/* V_eff (embedding) = V_xc (rho_total) + d T_nadd(rho_total, rho_subsys) / d rho_subsys
                                     = V_xc (rho_totoal) + d T (rho_total) / d rho_total - d T (rho_subsys) /d rho_subsys 
		 */

		/* XC-part */
		embedding_xc_functional->first (ds, factor, dp);

		/* nonadditive kinetic energy part */
		embedding_kin_functional->first (ds, factor, dp);
		embedding_kin_functional->first (ds, -factor, dp->subsystem[0]);
	}
}

static void
emb_second (FunSecondFuncDrv *ds, real factor, const FunDensProp *dp)
{
	quit_ ("Second functional derivatives for embedding are not implemented yet\n");
}

static void
emb_third (FunThirdFuncDrv *ds, real factor, const FunDensProp *dp)
{
	quit_ ("Third functional derivatives for embedding are not implemented yet\n");
}

Functional EmbeddingFunctional = {
	"embedding",
	emb_isgga,
	emb_read,
	emb_report,
	emb_energy,
	emb_first,
	emb_second,
	emb_third
};
 
/* 
 * EmbInteractionFunctional and enbint_* functions:
 *
 * This is a generic implementation of a "Embedding Interaction" functional,
 * as described above.
 *
 * It uses the same variables and initialization routine as EmbeddingFunctional,
 * so that the two can be interchanged in the middle of calculations without 
 * problems, or both can be used together (have thought about it in details yet,
 * maybe there will be a need to do that)
 *
 */

static void
embint_report (void)
{
	fort_print ("\n     Embedding Intercation Functional:");
	fort_print ("          Exchange-correlation functional: %20s", embedding_xc_functional->name);
	if (embedding_xc_functional->report) 
		embedding_xc_functional->report ();
	fort_print ("          Kinetic energy functional:       %20s \n", embedding_kin_functional->name);
	if (embedding_kin_functional->report)
		embedding_kin_functional->report ();
}

static real
embint_energy (const FunDensProp *dp)
{
	real energy = 0;

	if (!dp->subsystem[0] || !dp->subsystem[1]) {
		fort_print ("embint_first called but dp->subsystem_dp not initialized\n");
		fort_print ("This should not happen ... \n");
		quit_ ("Error in embint_first");
	} else {
		/* This is the XC- and kinetic energy part of the interaction energy between
		   regions I and II. All electrostatic interactions are NOT included           */

		/* E (int) = E_xc (int) + T_nadd
		           = E_xc (rho_tot) - E_xc (rho1) - E_xc (rho2) +
			     T (rho_tot) - T (rho1) - T (rho2)
		 */

		/* XC-part */
		energy += embedding_xc_functional->func (dp);
		energy -= embedding_xc_functional->func (dp->subsystem[0]);
		energy -= embedding_xc_functional->func (dp->subsystem[1]);

		/* nonadditive kinetic energy part */
		energy += embedding_kin_functional->func (dp);
		energy -= embedding_kin_functional->func (dp->subsystem[0]);
		energy -= embedding_kin_functional->func (dp->subsystem[1]);

	}
	return energy;
}

static void
embint_first (FunFirstFuncDrv *ds, real factor, const FunDensProp *dp)
{
	if (!dp->subsystem[0] || !dp->subsystem[1]) {
		fort_print ("embint_first called but dp->subsystem_dp not initialized\n");
		fort_print ("This should not happen ... \n");
		quit_ ("Error in embint_first");
	} else {
		/* V_eff (interaction) = d E_xc(int) (rho_tot, rho_1) / d rho_1 + d T_nadd(rho_total, rho_1) / d rho_1
                                       = d E_xc (rho_total) / d rho_total - d E_xc (rho_1) /d rho_1  
                                         + d T (rho_total) / d rho_total - d T (rho_1) /d rho_1  
		 */

		/* XC-part */
		embedding_xc_functional->first (ds, factor, dp);
		embedding_xc_functional->first (ds, -factor, dp->subsystem[0]);

		/* nonadditive kinetic energy part */
		embedding_kin_functional->first (ds, factor, dp);
		embedding_kin_functional->first (ds, -factor, dp->subsystem[0]);
	}
}

static void
embint_second (FunSecondFuncDrv *ds, real factor, const FunDensProp *dp)
{
	quit_ ("Second functional derivatives for embedding "
               "are not implemented yet\n");
}

static void
embint_third (FunThirdFuncDrv *ds, real factor, const FunDensProp *dp)
{
	quit_ ("Third functional derivatives for embedding "
               "are not implemented yet\n");
}

Functional EmbInteractionFunctional = {
	"emb_interaction",
	emb_isgga,
	emb_read,
	embint_report,
	embint_energy,
	embint_first,
	embint_second,
	embint_third
};
