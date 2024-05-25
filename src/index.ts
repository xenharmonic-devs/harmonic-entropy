import {gcd, valueToCents} from 'xen-dev-utils';
import {conv, padded64} from './helpers';

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

  // Size of the padded kernel
  const kernelSize = Math.ceil((max - min) / res) + 1;

  // build kernel (with additional padding to avoid circular convolution issues)
  const k = padded64(kernelSize);
  const ak = padded64(kernelSize);

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

    const icompl = 1 / rcompl;
    const acompl = Math.pow(rcompl, -a);

    //start building kernel, first check for rounded off case that doesn't need interpolation
    if (!mu) {
      k[index] += icompl;
      ak[index] += acompl;
    }
    //or else we do need interpolation
    else {
      k[index] += icompl * (1 - mu);
      k[index + 1] += icompl * mu;

      ak[index] += acompl * (1 - mu);
      ak[index + 1] += acompl * mu;
    }
  }

  // do convolution

  // now work out gaussian
  const g = padded64(kernelSize);
  const ag = padded64(kernelSize);
  let g_sum = 0;
  const s = -1 / (2 * scents * scents);
  for (let i = g.length - 1; i >= 0; i--) {
    const c = i * res + min;
    const gval =
      (1 / (scents * 2 * Math.PI)) * Math.exp(Math.pow(c - min, 2) * s) +
      (1 / (scents * 2 * Math.PI)) *
        Math.exp(Math.pow(c - (g.length * res + min), 2) * s);
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
      const max = Math.floor(Math.sqrt(n));
      for (let i = 1; i <= max; i++) {
        const m = n / i;
        if (Number.isInteger(m) && gcd(i, m) === 1) {
          r.push([i, m]); //does numerator on left, denominator on right
          if (m !== i) r.push([m, i]);
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
