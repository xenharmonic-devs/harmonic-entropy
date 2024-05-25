import {describe, bench, beforeAll} from 'vitest';
import {precalculateRatios, HarmonicEntropyOptions, harmonicEntropy} from '..';

// These benchmarks are not comparative.
// They're here to inform optimization of the code from version to version.

describe('Ratio pre-calculator', () => {
  bench('Tenney', () => {
    precalculateRatios({
      N: 10000,
      series: 'tenney',
    } as HarmonicEntropyOptions);
  });

  bench('Farey', () => {
    precalculateRatios({
      N: 1000,
      series: 'farey',
    } as HarmonicEntropyOptions);
  });
});

const options: HarmonicEntropyOptions = {
  N: 10000,
  a: 1,
  s: 0.01,
  series: 'tenney',
  minCents: 0,
  maxCents: 2400,
  res: 1,
  normalize: false,
};
let r: [number, number][];
beforeAll(() => {
  r = precalculateRatios(options);
});

describe('Harmonic entropy calculator', () => {
  bench('Tenney', () => {
    harmonicEntropy(options, r);
  });
});
