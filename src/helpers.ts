import {ceilPow2, fft, ifftReal} from 'frost-fft';

// Internal helpers not part of the public API

export function padded64(length: number) {
  const l = ceilPow2(length);
  if (l === length) {
    return new Float64Array(l * 2);
  }
  return new Float64Array(l);
}

export function conv(a: Float64Array, b: Float64Array) {
  const [f_a_real, f_a_imag] = fft(a);
  const [f_b_real, f_b_imag] = fft(b);

  // Normalize result
  const inorm = 1 / a.length;

  // (a+bi)(c+di)
  // (ac - bd) + (ad + bc)i

  // Do the multiplication out
  for (let i = a.length - 1; i >= 0; i--) {
    // Re-use arrays
    const temp = f_a_real[i];
    f_a_real[i] =
      (f_a_real[i] * f_b_real[i] - f_a_imag[i] * f_b_imag[i]) * inorm;
    f_a_imag[i] = (temp * f_b_imag[i] + f_a_imag[i] * f_b_real[i]) * inorm;
  }

  return ifftReal(f_a_real, f_a_imag);
}
