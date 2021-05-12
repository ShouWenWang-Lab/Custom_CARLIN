/*========================================================================
 * Aligns sequences according to CRISPR-Cas9 cutsites
 *
 * Can be used as a standalone file for other groups who may have 
 * their own CRISPR-Cas9 barcoding system.
 *
 * Author: Duluxan Sritharan. Hormoz Lab. Harvard Medical School.
 *
 * If you use this code, please cite:
 * 
 * S. Bowling, D. Sritharan, F. G. Osorio, M. Nguyen, P. Cheung, 
 * A. Rodiguez-Fraticelli, S. Patel, W-C. Yuan, Y. Fujiwara, B. E. Li, S. H. Orkin, 
 * S. Hormoz, F. D. Camargo. "An Engineered CRISPR-Cas9 Mouse Line for 
 * Simultaneous Readout of Lineage Histories and Gene Expression Profiles 
 * in Single Cells." Cell (2020), https://doi.org/10.1016/j.cell.2020.04.048
 *
 * Five Inputs:
 *
 * 1. Sequence to be aligned against reference in UINT8 format as produced 
 *    by nt2int: (A = 1, C = 2, G = 3, T = 4). Should have no gaps.
 *
 * 2. Reference template against which the observed sequence (arg 1) is to 
 *    be aligned, in the same format. Should have no gaps.
 *    
 * 3. Double array of opening gap penalties. Values should be low at sites
 *    where CRISPR-Cas9 is expected to cut, and high otherwise. See values
 *    used in CARLIN as an example.
 *    
 * 4. Double array of closing gap penalties. Values should be low at sites
 *    where CRISPR-Cas9 is expected to cut, and high otherwise. See values
 *    used in CARLIN as an example.
 *
 * 5. 4x4 double array of match/mismatch scores. CARLIN uses NUC44 values.
 *
 * Two Outputs:
 *
 * 1. Score of optimal alignment
 *
 * 2. 2xN UINT8 array of aligned sequence with top row corresponding to the 
 *    input sequence, and the bottom to the reference sequence, following
 *    the same format as the reference and input sequence. Gaps are marked 
 *    as 0s.
 *
 * This function was originally written to be called from a MATLAB wrapper, 
 * so there is limited input validation. You might trip on this if you use
 * this code as a standalone file:
 *
 * All input values are expected to be 1-indexed as in the manuscript, and 
 * should be padded to the left by 0. The exception is the opening penalty,
 * where, as in the manuscript, we reserve 0th index for the leading insertion
 * penalty. This ends up being convenient so that the input argument can
 * be used directly instead of reallocating space just to shift the
 * array. The output doesn't follow this convention, and strings are 0-indexed.
 * See cas_align_mex.cpp for a sample invocation.
 *
 *======================================================================*/

#ifndef __CAS9_ALIGN_HPP
#define __CAS9_ALIGN_HPP

#include<cstdlib>
#include<vector>

double cas9_align(const unsigned char *seq, const size_t Lseq, 
                  const unsigned char *ref, const size_t Lref,
                  const double* open_penalty, const double* close_penalty,
                  const double* sub_score, std::vector<unsigned char> &al_seq, std::vector<unsigned char> &al_ref);

#endif