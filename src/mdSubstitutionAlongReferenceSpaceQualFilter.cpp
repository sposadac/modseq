#include <Rcpp.h>

using namespace Rcpp;

double sumVec(NumericVector x, int i) {
  double total = 0;
  for(int k = 0; k <= i; ++k) {
    total += x[k];
  }
  return total;
}

// [[Rcpp::export]] 

CharacterVector mdSubstitutionAlongReferenceSpaceQualFilter(NumericVector md_ops_lengths, 
std::vector<std::string> md_ops, CharacterVector cigar_ops, NumericVector cigar_ops_lengths, 
std::string read, std::string qual, NumericVector mod_len_cumsum, CharacterVector names_mod_len_cumsum, 
int offset, int qthold = 15) {
  if (md_ops_lengths.size() == 1) {
    return CharacterVector::create();
  } else {
    std::string* tmp_ret;
    tmp_ret = new std::string [md_ops_lengths.size()];
    int countInsertion = 0;
    int cum_LenInsertion = 0;
    int* posInsertion;
    int* lenInsertion;
    int max = cigar_ops_lengths.size();
    posInsertion = new int [max];
    lenInsertion = new int [max];

    for (int k = 0; k < max; k++) {
      if (cigar_ops[k] == "I") {
        if (cigar_ops[0] == "S") {
          posInsertion[countInsertion] = cigar_ops_lengths[1] + 1;
          for (int a = 2; a < k ; a++) {
            posInsertion[countInsertion] = posInsertion[countInsertion] + cigar_ops_lengths[a];
          }
        } else {
          posInsertion[countInsertion] = cigar_ops_lengths[0] + 1;
          for (int a = 1; a < k ; a++) {
            posInsertion[countInsertion] = posInsertion[countInsertion] + cigar_ops_lengths[a];
          }
        }
        lenInsertion[countInsertion] = cigar_ops_lengths[k];
        if (countInsertion >= 1) {
            cum_LenInsertion = cum_LenInsertion + lenInsertion[countInsertion - 1];
            posInsertion[countInsertion] = posInsertion[countInsertion] - cum_LenInsertion;
        } 
        countInsertion++;
      }
    }
  
    int count = 0;
    int countPos = 0;
    int pos = 0;
    int aux_i = 0;
    int qualInt = 0;
    
    for (int i = 0; i < (md_ops_lengths.size()-1); i++) {
      if (md_ops[i+1].substr(0,1) != "^") {
        double aux_sum = sumVec(md_ops_lengths, i);
        double index = aux_sum + count; //w.r.t trimmed read
        if (countInsertion >= 1) {
          for (int j = 0; j < countInsertion; j++) {
            if (index >= posInsertion[j]) {
              index += lenInsertion[j];
            }
          }
        }
        qualInt = (int)qual[index] - 33;
        
        if (qualInt >= qthold) {
          std::string substitution;
          substitution = read[index];
          pos = aux_sum + countPos + offset; //w.r.t reference sequence
          int aux_mod = 0;
          for (int it = 0; it < mod_len_cumsum.size(); it++) {
            if (mod_len_cumsum[it] - 1 <= pos) {
              aux_mod ++;
            } else {
              break; //TODO: verify Break
            }
          }
          std::string mod;
          mod = names_mod_len_cumsum[aux_mod - 1];
          std::ostringstream aux_modPos;
          aux_modPos << (pos - mod_len_cumsum[aux_mod - 1] + 2);
          std::string modPos;
          modPos = aux_modPos.str();
          tmp_ret[aux_i] = mod + ":" + modPos + ":" + substitution;
          aux_i++;
        }
      } else if (md_ops[i+1].substr(0,1) == "^") {
        for (int j = 0; j < countInsertion; j++) {
          if ((sum(head(md_ops_lengths, i)) + count) <  posInsertion[j]) {
            posInsertion[j] = posInsertion[j] - md_ops[i+1].size() - 1; //'^' character should be discounted
          }
        }
        count --;
        countPos = countPos + md_ops[i+1].size() - 2; //one was added, by initialization, plus '^' character should be discounted
      }
      count ++;
      countPos ++;
    }

    CharacterVector ret(aux_i);
    for (int i = 0; i < aux_i; i++) {
      ret[i] = tmp_ret[i];
    }
    delete [] posInsertion;
    delete [] lenInsertion;
    delete [] tmp_ret;
    return ret;    
  }
}