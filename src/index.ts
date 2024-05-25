import {Fraction, FractionValue, gcd, valueToCents} from 'xen-dev-utils';
import {conv, padded64} from './helpers';

/**
 * Options for configuring harmonic entropy calculation.
 */
export type HarmonicEntropyOptions = {
  /** Max height of rationals (Benedetti (default 10k) or Wilson (default 1000) depending on series) */
  N?: number;
  /** Gaussian frequency deviation (default 0.01) */
  s?: number;
  /** Rényi order (default 1.0 i.e Shanon entropy) */
  a?: number;
  /** Series of rationals to use (default 'tenney') */
  series?: 'tenney' | 'farey';
  /** Lower bound of tabulation (default 0) */
  minCents?: number;
  /** Upper bound of tabulation (default 2400) */
  maxCents?: number;
  /** Tabulation delta in cents (default 1) */
  res?: number;
  /** Flag to normalize the result by Hartley entropy (default false) */
  normalize?: boolean;
};

/**
 * Calculate a harmonic entropy table.
 * @param options Parameters for calculating the table.
 * @param ratios Ratios precalculated using {@link precalculateRatios}.
 * @returns Array of [cents, entropy (in natural units)].
 */
export function harmonicEntropy(
  options: HarmonicEntropyOptions,
  ratios: number[][]
): [number, number][] {
  const scents = valueToCents((options.s ?? 0.01) + 1);
  let a = options.a ?? 1;
  const res = options.res ?? 1;
  const series = options.series ?? 'tenney';
  const normalize = !!options.normalize;
  let min = options.minCents ?? 0;
  let max = options.maxCents ?? 2400;

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

  for (let i = ratios.length - 1; i >= 0; i--) {
    let rcompl;

    const rcent = valueToCents(ratios[i][0] / ratios[i][1]);

    if (series === 'tenney') rcompl = Math.sqrt(ratios[i][0] * ratios[i][1]);
    else if (series === 'farey') rcompl = ratios[i][1];
    else throw new Error(`Unsupported series ${series}`);

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

/**
 * Precalculate the set of ratios considered when calculating {@link harmonicEntropy}.
 * @param options Parameters for calculating the table.
 * @returns Array of [numerator, denominator] pairs.
 */
export function precalculateRatios(options: HarmonicEntropyOptions) {
  const r = new Array<[number, number]>();

  let n = options.N;
  const series = options.series ?? 'tenney';
  if (series === 'tenney') {
    n ??= 10000;
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
  } else if (series === 'farey') {
    n ??= 1000;
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

/**
 * Construct a harmonic entropy calculator for individual musical intervals.
 */
export class EntropyCalculator {
  private options: HarmonicEntropyOptions;
  private ratios?: [number, number][];
  table: [number, number][];

  constructor(options?: HarmonicEntropyOptions) {
    this.options = {...options};
    // Hack to disable calculation during revification
    const series = this.options.series;
    if (series === undefined || series === 'tenney' || series === 'farey') {
      this.ratios = precalculateRatios(this.options);
      this.table = harmonicEntropy(this.options, this.ratios);
    } else {
      // Keep TypeScript happy.
      this.table = [];
    }
  }

  /**
   * Revive a {@link EntropyCalculator} instance produced by `EntropyCalculator.toJSON()`. Return everything else as is.
   *
   * Intended usage:
   * ```ts
   * const data = JSON.parse(serializedData, EntropyCalculator.reviver);
   * ```
   *
   * @param key Property name.
   * @param value Property value.
   * @returns Deserialized {@link EntropyCalculator} instance or other data without modifications.
   */
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  static reviver(key: string, value: any) {
    if (
      typeof value === 'object' &&
      value !== null &&
      value.type === 'EntropyCalculator'
    ) {
      const series = value.options.series;
      const options = {...value.options};
      options.series = '__revived';
      const instance = new EntropyCalculator(options);
      // Bypass setter to prevent recalculation
      instance.options.series = series;
      const minCents = instance.minCents;
      const res = instance.res;
      instance.table = value.tableY.map((y: number, i: number) => [
        minCents + i * res,
        y,
      ]);
      return instance;
    }
    return value;
  }

  /**
   * Serialize the entropy calculator to a JSON compatible object.
   * @returns The serialized object with property `type` set to `'EntropyCalculator'`.
   */
  toJSON() {
    return {
      type: 'EntropyCalculator',
      options: this.options,
      tableY: this.table.map(xy => xy[1]),
    };
  }

  private recalculate() {
    if (!this.ratios) {
      this.ratios = precalculateRatios(this.options);
    }
    this.table = harmonicEntropy(this.options, this.ratios);
  }

  get N() {
    if (this.series === 'tenney') {
      return 10000;
    }
    return 1000;
  }
  set N(value: number) {
    this.options.N = value;
    this.ratios = undefined;
    this.recalculate();
  }

  get a() {
    return this.options.a ?? 1.0;
  }
  set a(value: number) {
    this.options.a = value;
    this.recalculate();
  }

  get s() {
    return this.options.s ?? 0.01;
  }
  set s(value: number) {
    this.options.s = value;
    this.recalculate();
  }

  get series() {
    return this.options.series ?? 'tenney';
  }
  set series(value: 'tenney' | 'farey') {
    this.options.series = value;
    this.recalculate();
  }

  get minCents() {
    return this.options.minCents ?? 0;
  }
  set minCents(value: number) {
    this.options.minCents = value;
    this.recalculate();
  }

  get maxCents() {
    return this.options.maxCents ?? 2400;
  }
  set maxCents(value: number) {
    this.options.maxCents = value;
    this.recalculate();
  }

  get res() {
    return this.options.res ?? 1;
  }
  set res(value: number) {
    this.options.res = value;
    this.recalculate();
  }

  get normalize() {
    return !!this.options.normalize;
  }
  set normalize(value: boolean) {
    this.options.normalize = value;
    this.recalculate();
  }

  /**
   * Calculate the harmonic entropy of a rational number.
   * @param value Fractional value.
   * @returns The harmonic entropy of the input in natural units.
   */
  ofFraction(value: FractionValue) {
    let cents: number;
    if (typeof value === 'number') {
      cents = valueToCents(value);
    } else {
      cents = valueToCents(new Fraction(value).valueOf());
    }
    return this.ofCents(cents);
  }

  /**
   * Calculate the harmonic entropy of a musical interval measured in cents.
   * @param cents Width of the interval.
   * @returns The harmonic entropy of the input in natural units.
   */
  ofCents(cents: number) {
    if (isNaN(cents)) {
      throw new Error('Invalid input');
    }
    // Dyadic entropy is symmetric
    cents = Math.abs(cents);
    if (cents < this.minCents || cents > this.maxCents) {
      throw new Error('Value out of tabulated range');
    }
    let mu = (cents - this.minCents) / this.res;
    const index = Math.floor(mu);
    mu -= index;
    if (!mu) {
      return this.table[index][1];
    }
    return this.table[index][1] * (1 - mu) + this.table[index + 1][1] * mu;
  }
}
