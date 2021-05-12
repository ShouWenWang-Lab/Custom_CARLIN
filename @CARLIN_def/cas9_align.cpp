#include"cas9_align.hpp"
#include<limits>
#include<cassert>
#include<algorithm>
   
enum Mutation { S, D, I, N_muts };

inline double max3(const double* vals, bool* argmax3) {
     const double maxval = *std::max_element(vals, vals+3);
     argmax3[S] = (vals[S] == maxval);
     argmax3[D] = (vals[D] == maxval);
     argmax3[I] = (vals[I] == maxval);
     return maxval;
}

double cas9_align(const unsigned char *seq, const size_t Lseq, 
                  const unsigned char *ref, const size_t Lref,
                  const double* open_penalty, const double* close_penalty,
                  const double* sub_score, std::vector<unsigned char> &al_seq, std::vector<unsigned char> &al_ref)
{   
    // Subalignments ending with a substitution, deletion and insertion
    std::vector<std::vector<std::vector<double> > > score(N_muts, std::vector<std::vector<double> >(\
                                                          Lseq+1, std::vector<double>(\
                                                          Lref+1, -std::numeric_limits<double>::infinity())));
    
    // Precursor to optimal substitution, deletion and insertion subalignment
    std::vector<std::vector<std::vector<std::vector<bool> > > > backtrack(N_muts, std::vector<std::vector<std::vector<bool> > >(\
                                                                          N_muts, std::vector<std::vector<bool> >(\
                                                                          Lseq+1, std::vector<bool>(\
                                                                          Lref+1, false))));
    
    score[S][0][0] = 0;
    for (int j = 1; j <= Lseq; j++) {
        score[I][j][0] = -open_penalty[0];
        backtrack[I][I][j][0] = 1;
    }
    for (int k = 1; k <= Lref; k++) {
        score[D][0][k] = -open_penalty[0];
        backtrack[D][D][0][k] = 1;
    }   
    
    for (int j = 1; j <= Lseq; j++) {
        for (int k = 1; k <= Lref; k++) {
            
            // S: Create a subalignment(j,k) ending with a substitution by:
            //  1) Extending a substitution
            //  2) Ending a deletion
            //  3) Ending an insertion
            
            // D: Create a subalignment(j,k) ending with a deletion by:
            //  1) Opening up a deletion
            //  2) Adding a nucleotide on the reference only, extending an ongoing deletion
            //  3) Ending an insertion, and beginning a deletion (nop if del_penalty(k) > mismatch(j,k))
             
            // I: Create a subalignment(j,k) ending with an insertion by:
            //   1) Opening up an insertion
            //   2) Ending a deletion, and beginning an insertion (nop if del_penalty(k) > mismatch(j,k))
            //   3) Adding a nucleotide on the sequence only, extending an ongoing insertion
            
            double ss = sub_score[seq[j]*5+ref[k]];
            double score_jk[N_muts][N_muts] = \
                    {{score[S][j-1][k-1]+ss             , score[D][j-1][k-1]-close_penalty[k-1]+ss, score[I][j-1][k-1]-close_penalty[k]+ss},
                     {score[S][j  ][k-1]-open_penalty[k], score[D][j  ][k-1]                      , score[I][j  ][k-1]                      },
                     {score[S][j-1][k  ]-open_penalty[k], score[D][j-1][k  ]                      , score[I][j-1][k  ]                      }};
             
             for (int m = 0; m < N_muts; m++) {
                 bool out[N_muts];
                 score[m][j][k] = max3(score_jk[m], out);
                 backtrack[S][m][j][k] = out[S];
                 backtrack[D][m][j][k] = out[D];
                 backtrack[I][m][j][k] = out[I];
             }
        }
    }
    
    double best_score;
    int cur_state;
    
    {
        double score_LL[N_muts] = {score[S][Lseq][Lref], score[D][Lseq][Lref], score[I][Lseq][Lref]};
        bool state_LL[N_muts];
        best_score = max3(score_LL, state_LL);
        cur_state = (state_LL[S] ? S : (state_LL[D] ? D : (state_LL[I] ? I : -1)));
    }
    
    int cur_j = Lseq;
    int cur_k = Lref;
    
    while (cur_j || cur_k) {
        
        assert(cur_state >= 0 && cur_state < N_muts);
        
        const bool cur_bt[N_muts] = {backtrack[S][cur_state][cur_j][cur_k],backtrack[D][cur_state][cur_j][cur_k],backtrack[I][cur_state][cur_j][cur_k]};        
        
        if (cur_state == S){
            al_seq.push_back(seq[cur_j]);
            al_ref.push_back(ref[cur_k]);
            cur_j = cur_j-1;
            cur_k = cur_k-1;
            cur_state = ( cur_bt[S] ? S : ( cur_bt[D] ? D : ( cur_bt[I] ? I : -1 ) ) );                
        }
        else if(cur_state == D) {
            al_seq.push_back(0);
            al_ref.push_back(ref[cur_k]);
            cur_k = cur_k-1;
            cur_state = ( cur_bt[D] ? D : ( cur_bt[S] ? S : ( cur_bt[I] ? I : -1 ) ) );                
        }
        else if (cur_state == I) {
            al_seq.push_back(seq[cur_j]);
            al_ref.push_back(0);
            cur_j = cur_j-1;
            cur_state = ( cur_bt[I] ? I : ( cur_bt[S] ? S : ( cur_bt[D] ? D : -1 ) ) );
        }
        else {
            cur_state = -1;        
        }
    }
    
    std::reverse(al_seq.begin(), al_seq.end());
    std::reverse(al_ref.begin(), al_ref.end());
    
    return best_score;
}
