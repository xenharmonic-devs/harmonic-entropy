import {describe, it, expect} from 'vitest';
import {
  HarmonicEntropyInfo,
  conv,
  harmonicEntropy,
  preCalcRatios,
} from '../index';

describe('Convolution', () => {
  it('causes no change to a single impulse', () => {
    const impulse = [1, 0, 0, 0, 0];
    const result = conv(impulse, impulse);
    expect(result).toEqual(new Float64Array([1, 0, 0, 0, 0]));
  });

  it('convolvest two arrays', () => {
    const a = [0.0, 1.0, 1.0, 0, 0, 0, 0, 0];
    const b = [0.5, 1.0, 0.5, 0, 0, 0, 0, 0];
    const result = conv(a, b);
    expect(result.map(n => Math.round(n * 1024) / 1024)).toEqual(
      new Float64Array([0, 0.5, 1.5, 1.5, 0.5, 0, 0, 0])
    );
  });
});

describe('Ratio pre-calculator', () => {
  it('calculates ratios for tenney series', () => {
    const r = preCalcRatios({
      N: 10,
      series: 'tenney',
    } as HarmonicEntropyInfo);
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
    const r = preCalcRatios({
      N: 5,
      series: 'farey',
    } as HarmonicEntropyInfo);
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

describe('Harmonic entropy calculator', () => {
  it('calculates harmonic entropy', () => {
    const HEinfo: HarmonicEntropyInfo = {
      N: 1000,
      a: 1,
      s: 0.01,
      series: 'tenney',
      dist: 'log',
      mincents: 0,
      maxcents: 1200,
      res: 1,
      normalize: false,
    };
    const r = preCalcRatios(HEinfo);
    const entropy = harmonicEntropy(HEinfo, r);
    expect(entropy).toHaveLength(1201);
    let x: number, y: number;
    [x, y] = entropy[100];
    expect(x).toBe(100);
    expect(y).toBeCloseTo(2.755);
    [x, y] = entropy[700];
    expect(x).toBe(700);
    expect(y).toBeCloseTo(1.15239);
  });
});
