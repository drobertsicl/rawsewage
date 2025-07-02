#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
List rcpp_process_batch(RawVector raw_data,
                          NumericVector adds,
                          IntegerVector indices,
                          double batch_start,
                          Nullable<NumericVector> calib_coefs = R_NilValue,
                          bool use_integer_intensity = false,
                          bool use_float_mz = false) {
  
  List results(indices.size());
  NumericVector coefs;
  bool use_calib = false;
  if (calib_coefs.isNotNull()) {
    coefs = calib_coefs.get();
    use_calib = true;
  }
  
  for (int i = 0; i < indices.size(); i++) {
    int idx = indices[i] - 1;
    int start = adds[idx] - batch_start;
    int end = (idx + 1 < adds.size()) ? adds[idx + 1] - batch_start : raw_data.size();
    if (start < 0 || end > raw_data.size()) continue;
    
    RawVector byte = raw_data[Range(start, end - 1)];
    int n_peaks = byte.size() / 8;

    // Allocate vectors
    std::vector<float> mz_vec(n_peaks);
    std::vector<double> mz_dvec(n_peaks);
    std::vector<int> intensity_ivec(n_peaks);
    std::vector<double> intensity_dvec(n_peaks);

    for (int j = 0; j < n_peaks; ++j) {
      int pos = j * 8;
      uint8_t b1 = byte[pos];
      uint8_t b2 = byte[pos + 1];
      uint8_t b3 = byte[pos + 2];
      uint8_t b4 = byte[pos + 3];
      uint8_t micro = byte[pos + 4];
      uint8_t frac  = byte[pos + 5];
      uint8_t add   = byte[pos + 6];
      uint8_t mz_b  = byte[pos + 7];
      
      uint8_t masked = mz_b & 0x7F;
      int num4bit = masked >> 3;
      int num3bit = masked & 0x07;
      
      double base_mz = 8 * std::pow(2, num4bit - 6) * num3bit;
      double add_mz = add * std::pow(2, num4bit - 11);
      double frac_mz = frac / std::pow(2, 19 - num4bit);
      double micro_mz = micro / std::pow(2, 27 - num4bit);
      
      double mz_val = base_mz + add_mz + frac_mz + micro_mz;

      // Calibration
      if (use_calib) {
        double root_mz = std::sqrt(mz_val);
        double acc = 1.0, total = 0.0;
        for (int c = 0; c < coefs.size(); ++c) {
          total += coefs[c] * acc;
          acc *= root_mz;
        }
        mz_val = total * total;
      }

      if (use_float_mz) {
        mz_vec[j] = static_cast<float>(mz_val);
      } else {
        mz_dvec[j] = mz_val;
      }

      // Intensity calculation
      double frac1 = 0, frac2 = 0, frac3 = 0;
      for (int k = 0; k < 8; ++k) {
        if (b1 & (1 << (7 - k))) frac1 += std::pow(2, -(14 + k));
        if (b2 & (1 << (7 - k))) frac2 += std::pow(2, -(6 + k));
      }
      int mult = 1;
      if (b3 & (1 << 7)) mult *= 4;
      if (b3 & (1 << 6)) mult *= 2;
      bool subtract1 = b3 & (1 << 5);
      
      for (int bit_pos = 4; bit_pos <= 8; ++bit_pos) {
        if (b3 & (1 << (8 - bit_pos))) frac3 += std::pow(2, -(bit_pos - 3));
      }

      double temp_sum = frac1 + frac2 + frac3;
      if (subtract1) temp_sum -= 1.0;

      int exp_factor = 0;
      if (b4 & 8) exp_factor += 32;
      if (b4 & 4) exp_factor += 16;
      if (b4 & 2) exp_factor += 8;
      if (b4 & 1) exp_factor += 4;

      double intensity_val = temp_sum * mult * std::pow(2.0, exp_factor);

      if (use_integer_intensity) {
        intensity_ivec[j] = static_cast<int>(std::round(intensity_val));
      } else {
        intensity_dvec[j] = intensity_val;
      }
    }

    // Create output
    List entry;
    if (use_float_mz) {
      NumericVector mz_out(n_peaks);
      for (int i = 0; i < n_peaks; i++) mz_out[i] = static_cast<double>(mz_vec[i]);
      entry["mz"] = mz_out;
    } else {
      entry["mz"] = wrap(mz_dvec);
    }

    if (use_integer_intensity) {
      entry["intensity"] = wrap(intensity_ivec);
    } else {
      entry["intensity"] = wrap(intensity_dvec);
    }

    results[i] = entry;
  }

  return results;
}
