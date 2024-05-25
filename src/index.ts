import {fft, ifftReal} from 'frost-fft';
import {gcd, valueToCents} from 'xen-dev-utils';

export type HarmonicEntropyInfo = {
  N: number;
  s: number;
  a: number;
  series: 'tenney' | 'farey';
  dist: 'linear' | 'log';
  mincents: number;
  maxcents: number;
  res: number;
  normalize: boolean;
};

export function conv(olda: number[], oldb: number[]) {
  if (olda.length !== oldb.length) {
    throw new Error('conv(...): input arrays must be the same length');
  }

  const a = olda.concat();
  const b = oldb.concat();

  const len = a.length;
  let minlen = 1;
  while (minlen < len) {
    // 2 is to avoid circular convolution distortion
    minlen *= 2;
  }

  for (let i = len; i < minlen; i++) {
    a[i] = 0;
    b[i] = 0;
  }

  const a_float64 = new Float64Array(a);
  const b_float64 = new Float64Array(b);
  const [f_a_real, f_a_imag] = fft(a_float64);
  const [f_b_real, f_b_imag] = fft(b_float64);

  // Reuse arrays
  const f_out_real = a_float64;
  const f_out_imag = b_float64;
  // Normalize result
  const inorm = 1 / minlen;

  // (a+bi)(c+di)
  // (ac - bd) + (ad + bc)i

  // Do the multiplication out
  for (let i = minlen - 1; i >= 0; i--) {
    f_out_real[i] =
      (f_a_real[i] * f_b_real[i] - f_a_imag[i] * f_b_imag[i]) * inorm;
    f_out_imag[i] =
      (f_a_real[i] * f_b_imag[i] + f_a_imag[i] * f_b_real[i]) * inorm;
  }

  return ifftReal(f_out_real, f_out_imag).slice(0, len);
}

export function harmonicEntropy(HEinfo: HarmonicEntropyInfo, locr: number[][]) {
  //globals
  const {res, dist, series, normalize} = HEinfo;
  let a = HEinfo.a;
  const scents = valueToCents(HEinfo.s + 1);
  let min = HEinfo.mincents;
  let max = HEinfo.maxcents;

  if (dist === 'linear') {
    min = 0;
    max = 1200;
  }

  const outLen = Math.ceil((max - min) / res) + 1;

  const arrayPadding = Math.ceil((100 * scents) / res);
  const padding = arrayPadding * res;
  min -= padding;
  max += padding;

  a = a === 1 ? (a = 1.0000000001) : a;
  let rcount = 0;

  // build kernel
  const k = new Array<number>();
  const ak = new Array<number>();

  for (let i = Math.ceil((max - min) / res); i >= 0; i--) {
    k[i] = 0;
    ak[i] = 0;
  }

  for (let i = locr.length - 1; i >= 0; i--) {
    let rcent, rcompl;

    if (dist === 'log') rcent = valueToCents(locr[i][0] / locr[i][1]);
    else if (dist === 'linear') rcent = (1200 * locr[i][0]) / locr[i][1];
    else throw new Error(`Unsupported dist ${dist}`);

    if (series === 'tenney') rcompl = Math.sqrt(locr[i][0] * locr[i][1]);
    else if (series === 'farey') rcompl = locr[i][1];
    else throw new Error(`Unsupported seriese ${series}`);

    //check for bounds to optimize
    if (rcent < min || rcent > max) continue;

    rcount++;

    let mu = (rcent - min) / res;
    const index = Math.floor(mu);
    mu -= index;

    //start building kernel, first check for rounded off case that doesn't need interpolation
    if (!mu) {
      k[index] += 1 / rcompl;
      ak[index] += 1 / Math.pow(rcompl, a);
    }
    //or else we do need interpolation
    else {
      k[index] += (1 / rcompl) * (1 - mu);
      k[index + 1] += (1 / rcompl) * mu;

      ak[index] += (1 / Math.pow(rcompl, a)) * (1 - mu);
      ak[index + 1] += (1 / Math.pow(rcompl, a)) * mu;
    }
  }

  // do convolution
  // first pad to a power of two to make the convolution easier for numerous reasons
  let minlen = 1;
  while (minlen < 2 * k.length)
    // 2 is to avoid circular convolution distortion
    minlen *= 2;

  for (let i = k.length; i < minlen; i++) {
    k[i] = 0;
    ak[i] = 0;
  }

  // now work out gaussian
  const g = new Array<number>(minlen);
  const ag = new Array<number>(minlen);
  let g_sum = 0;
  for (let i = minlen - 1; i >= 0; i--) {
    const c = i * res + min;
    const gval =
      (1 / (scents * 2 * Math.PI)) *
        Math.exp(-(Math.pow(c - min, 2) / (2 * scents * scents))) +
      (1 / (scents * 2 * Math.PI)) *
        Math.exp(
          -(Math.pow(c - (minlen * res + min), 2) / (2 * scents * scents))
        );
    g[i] = gval;
    g_sum += gval;
  }
  // normalize gaussian
  for (let i = g.length - 1; i >= 0; i--) {
    g[i] /= g_sum;
    ag[i] = Math.pow(g[i], a);
  }

  const ent = conv(ak, ag);
  const nrm = conv(k, g);

  let nrmfct;
  if (normalize) {
    nrmfct = Math.log(rcount);
  } else nrmfct = 1;

  // trim answer and out
  const out = new Array<[number, number]>();
  for (let i = outLen - 1; i >= 0; i--) {
    const j = i + arrayPadding;
    const outval =
      ((1 / (1 - a)) * Math.log(ent[j] / Math.pow(nrm[j], a))) / nrmfct;
    out[i] = [j * res + min, outval];
  }
  return out;
}

export function preCalcRatios(HEinfo: HarmonicEntropyInfo) {
  const r = new Array<[number, number]>();

  let n = HEinfo.N;
  if (HEinfo.series === 'tenney') {
    do {
      for (let i = 1; i <= Math.floor(Math.sqrt(n)); i++) {
        if (n / i === Math.round(n / i) && gcd(i, n / i) === 1) {
          r.push([i, n / i]); //does numerator on left, denominator on right
          if (n / i !== i) r.push([n / i, i]);
        }
      }
    } while (--n >= 0);
  } else if (HEinfo.series === 'farey') {
    do {
      for (let i = 0; i <= n; i++) {
        if (gcd(i, n) === 1) {
          r.push([n, i]);
          if (n !== i) r.push([i, n]);
        }
      }
    } while (--n >= 0);
  }
  return r;
}
