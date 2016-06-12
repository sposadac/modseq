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

CharacterVector mdSubstitutionAlongReferenceSpace(NumericVector md_ops_lengths, 
std::vector<std::string> md_ops, CharacterVector cigar_ops, 
NumericVector cigar_ops_lengths, std::string read, NumericVector mod_len_cumsum,
CharacterVector names_mod_len_cumsum, int offset) {
  //CharacterVector ret = CharacterVector::create("0");
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
//          if(lenInsertion[countInsertion - 1] > 1){
            cum_LenInsertion = cum_LenInsertion + lenInsertion[countInsertion - 1];
            posInsertion[countInsertion] = posInsertion[countInsertion] - cum_LenInsertion;
//          }
        } 
        countInsertion++;
      }
	  }
    
    //ret[0] = (char) countInsertion; //check point
    int count = 0;
    int countPos = 0;
    int pos = 0;
    int aux_i = 0;
    
    for (int i = 0; i < (md_ops_lengths.size()-1); i++) {
      if (md_ops[i+1].substr(0,1) != "^") {
        // sugar function 'head' -- previous implementation
        /*
        int index = sum(head(md_ops_lengths, i+1)) + count; //w.r.t trimmed read
        int aux_sum = 0;
        for(int k = 0; k <= i; k++){
          aux_sum += md_ops_lengths[k];
        }*/
        double aux_sum = sumVec(md_ops_lengths, i);
        double index = aux_sum + count; //w.r.t trimmed read
        if (countInsertion >= 1) {
          for (int j = 0; j < countInsertion; j++) {
            if (index >= posInsertion[j]) {
              index += lenInsertion[j];
            }
          }
        }
        std::string substitution;
        //substitution = as<std::string>(read[index]);
        substitution = read[index];
        pos = aux_sum + countPos + offset; //w.r.t reference sequence
        //ret[i] = (char) pos; //check point
        int aux_mod = 0;
        for (int it = 0; it < mod_len_cumsum.size(); it++){
          if (mod_len_cumsum[it] - 1 <= pos) {
            aux_mod ++;
          } else {
            break; //TODO: Break loop
          }
        }
        std::string mod;
        //if(aux_mod == 0) aux_mod = 1;
        mod = names_mod_len_cumsum[aux_mod - 1];
        std::ostringstream aux_modPos;
        aux_modPos << (pos - mod_len_cumsum[aux_mod - 1] + 2);
        std::string modPos;
        modPos = aux_modPos.str();
        /*CharacterVector mod; 
        mod = names_mod_cumsumLen[aux_mod - 1];
        CharacterVector modPos;
        modPos = pos - mod_cumsumLen[aux_mod - 1] + 2; */
        tmp_ret[aux_i] = mod + ":" + modPos + ":" + substitution;
        aux_i++;
      } else if (md_ops[i+1].substr(0,1) == "^") {
        //ret = md_ops[i+1].substr(0,1); // check
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
    for (int i=0; i < aux_i; i++) {
      ret[i] = tmp_ret[i];
    }
    delete [] posInsertion;
    delete [] lenInsertion;
    delete [] tmp_ret;
    return ret;
  }
}


