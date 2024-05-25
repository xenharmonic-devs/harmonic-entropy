import {describe, it, expect} from 'vitest';
import {padded64, conv} from '../helpers';

describe('Convolution', () => {
  it('causes no change to a single impulse', () => {
    const N = 5;
    const impulse = padded64(N);
    impulse[0] = 1;
    const result = conv(impulse, impulse);
    expect(result.slice(0, N)).toEqual(new Float64Array([1, 0, 0, 0, 0]));
  });

  it('convolves two arrays', () => {
    const N = 8;
    const a = padded64(N);
    const b = padded64(N);
    a[1] = 1;
    a[2] = 1;
    b[0] = 0.5;
    b[1] = 1;
    b[2] = 0.5;
    const result = conv(a, b);
    expect(result.slice(0, N).map(n => Math.round(n * 1024) / 1024)).toEqual(
      new Float64Array([0, 0.5, 1.5, 1.5, 0.5, 0, -0, 0])
    );
  });
});
