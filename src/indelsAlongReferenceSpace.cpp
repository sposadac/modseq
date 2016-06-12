#include <Rcpp.h>
using namespace Rcpp;

double sumVec(NumericVector x, int i, int s) {
  double total = 0;
  for(int k = s; k < i; ++k) {
    total += x[k];
  }
  return total;
}

// [[Rcpp::export]]
CharacterVector indelsAlongReferenceSpace(CharacterVector cigar_ops, NumericVector cigar_ops_lengths, 
NumericVector mod_len_cumsum, CharacterVector names_mod_len_cumsum, int offset) {
  std::string* tmp_ret;
  tmp_ret = new std::string [cigar_ops_lengths.size()];
  int i = 0;
  int s = 0;
  int count = 0;
  int countIn = 0;
  if (cigar_ops[0] == "S") {
    s = 1;
    ++i;
  }
  for (int i; i < cigar_ops_lengths.size(); i++) {
    if (cigar_ops[i] == "I" | cigar_ops[i] == "D") {
      int len = cigar_ops_lengths[i];
      double pos;
      pos = sumVec(cigar_ops_lengths, i, s) - countIn + offset + 1;
      int aux_mod = 0;
      std::string type;
      std::ostringstream aux_modPos;
      if (cigar_ops[i] == "I") { 
        type = "In";
        countIn += len;
        for (int it = 0; it < mod_len_cumsum.size(); it++) {
          if (mod_len_cumsum[it] < pos) {
            aux_mod ++;
          } else {
            break;
          }
        }
        aux_modPos << (pos - mod_len_cumsum[aux_mod - 1]);
      } else if (cigar_ops[i] == "D") {
        type = "Del";
        for (int it = 0; it < mod_len_cumsum.size(); it++) {
          if (mod_len_cumsum[it] <= pos) {
            aux_mod ++;
          } else {
            break; 
          }
        }
        aux_modPos << (pos - mod_len_cumsum[aux_mod - 1] + 1);
      }
      std::string mod;
      mod = names_mod_len_cumsum[aux_mod - 1];

      std::string modPos;
      modPos = aux_modPos.str();
      std::ostringstream tmp;
      tmp << (len);
      std::string aux_len;
      aux_len = tmp.str();
      tmp_ret[count] = type + ":" + mod + ":" + modPos + ":" + aux_len;
      ++count;
    }
  }
  
  CharacterVector ret(count);
  for (int i = 0; i < count; i++) {
    ret[i] = tmp_ret[i];
  }
  delete [] tmp_ret;
  return ret;
}