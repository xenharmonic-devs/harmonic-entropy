import {describe, bench, beforeAll} from 'vitest';
import {preCalcRatios, HarmonicEntropyInfo, harmonicEntropy} from '..';

// These benchmarks are not comparative.
// They're here to inform optimization of the code from version to version.

describe('Ratio pre-calculator', () => {
  bench('Tenney', () => {
    preCalcRatios({
      N: 10000,
      series: 'tenney',
    } as HarmonicEntropyInfo);
  });

  bench('Farey', () => {
    preCalcRatios({
      N: 1000,
      series: 'farey',
    } as HarmonicEntropyInfo);
  });
});

const HEinfo: HarmonicEntropyInfo = {
  N: 10000,
  a: 1,
  s: 0.01,
  series: 'tenney',
  dist: 'log',
  mincents: 0,
  maxcents: 2400,
  res: 1,
  normalize: false,
};
let r: [number, number][];
beforeAll(() => {
  r = preCalcRatios(HEinfo);
});

describe('Harmonic entropy calculator', () => {
  bench('Tenney', () => {
    harmonicEntropy(HEinfo, r);
  });
});
