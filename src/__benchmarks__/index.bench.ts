import {describe, bench} from 'vitest';
import {preCalcRatios, HarmonicEntropyInfo} from '..';

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
