import {describe, it, expect} from 'vitest';
import {
  EntropyCalculator,
  HarmonicEntropyOptions,
  harmonicEntropy,
  precalculateRatios,
} from '../index';

describe('Ratio pre-calculator', () => {
  it('calculates ratios for tenney series', () => {
    const r = precalculateRatios({
      N: 10,
      series: 'tenney',
    } as HarmonicEntropyOptions);
    expect(r).toEqual([
      [1, 10],
      [10, 1],
      [2, 5],
      [5, 2],
      [1, 9],
      [9, 1],
      [1, 8],
      [8, 1],
      [1, 7],
      [7, 1],
      [1, 6],
      [6, 1],
      [2, 3],
      [3, 2],
      [1, 5],
      [5, 1],
      [1, 4],
      [4, 1],
      [1, 3],
      [3, 1],
      [1, 2],
      [2, 1],
      [1, 1],
    ]);
  });

  it('calculates ratios for farey series', () => {
    const r = precalculateRatios({
      N: 5,
      series: 'farey',
    } as HarmonicEntropyOptions);
    expect(r).toEqual([
      [5, 1],
      [1, 5],
      [5, 2],
      [2, 5],
      [5, 3],
      [3, 5],
      [5, 4],
      [4, 5],
      [4, 1],
      [1, 4],
      [4, 3],
      [3, 4],
      [3, 1],
      [1, 3],
      [3, 2],
      [2, 3],
      [2, 1],
      [1, 2],
      [1, 0], // ?!
      [0, 1],
      [1, 1],
    ]);
  });
});

describe('Harmonic entropy calculator (function)', () => {
  it('calculates harmonic entropy', () => {
    const options: HarmonicEntropyOptions = {
      N: 1000,
      a: 1,
      s: 0.01,
      series: 'tenney',
      minCents: 0,
      maxCents: 1200,
      res: 1,
      normalize: false,
    };
    const r = precalculateRatios(options);
    const entropy = harmonicEntropy(options, r);
    expect(entropy).toHaveLength(1201);
    let [x, y] = entropy[100];
    expect(x).toBe(100);
    expect(y).toBeCloseTo(2.755);
    [x, y] = entropy[700];
    expect(x).toBe(700);
    expect(y).toBeCloseTo(1.15239);
  });

  it('supports fractional resolution', () => {
    const options: HarmonicEntropyOptions = {
      N: 1000,
      a: 1,
      s: 0.01,
      series: 'tenney',
      minCents: 0,
      maxCents: 100,
      res: 0.5,
      normalize: false,
    };
    const r = precalculateRatios(options);
    const entropy = harmonicEntropy(options, r);
    let [x, y] = entropy[51];
    expect(x).toBe(25.5);
    expect(y).toBeCloseTo(0.50917);
    [x, y] = entropy[200];
    expect(x).toBe(100);
    expect(y).toBeCloseTo(2.755);
  });

  it('supports coarce resolution', () => {
    const options: HarmonicEntropyOptions = {
      N: 1000,
      a: 1,
      s: 0.01,
      series: 'tenney',
      minCents: 0,
      maxCents: 100,
      res: 2,
      normalize: false,
    };
    const r = precalculateRatios(options);
    const entropy = harmonicEntropy(options, r);
    const [x, y] = entropy[50];
    expect(x).toBe(100);
    expect(y).toBeCloseTo(2.755);
  });
});

describe('Harmonic entropy calculator (class)', () => {
  it('calculates the harmonic entropy of 700 cents and 100 cents', () => {
    const options: HarmonicEntropyOptions = {
      N: 1000,
      series: 'tenney',
      minCents: 0,
      maxCents: 1200,
    };
    const entropy = new EntropyCalculator(options);
    expect(entropy.ofCents(700)).toBeCloseTo(1.15239);
    expect(entropy.ofCents(100)).toBeCloseTo(2.755);
  });

  it('calculates the harmonic entropy of 3/2 and 9/8', () => {
    const options: HarmonicEntropyOptions = {
      N: 1000,
      series: 'tenney',
      minCents: 0,
      maxCents: 1200,
    };
    const entropy = new EntropyCalculator(options);
    expect(entropy.ofFraction('3/2')).toBeCloseTo(1.14437);
    expect(entropy.ofFraction(9 / 8)).toBeCloseTo(2.4649);
  });

  it('has default values', () => {
    const options: HarmonicEntropyOptions = {
      N: 1000,
    };
    const entropy = new EntropyCalculator(options);
    expect(entropy.a).toBe(1);
    expect(entropy.s).toBe(0.01);
    expect(entropy.series).toBe('tenney');
    expect(entropy.minCents).toBe(0);
    expect(entropy.maxCents).toBe(2400);
    expect(entropy.normalize).toBe(false);
  });
});
