# harmonic-entropy
Compute harmonic entropy of a musical interval

## Installation
```bash
npm i harmonic-entropy
```

## Documentation
Documentation is hosted at the project [Github pages](https://xenharmonic-devs.github.io/harmonic-entropy).

To generate documentation locally run:
```bash
npm run doc
```
### Other resources
You can read about how harmonic entropy relates to musical concodrance in the [XenWiki article](https://en.xen.wiki/w/Harmonic_entropy).

## Examples
Compute the entropy graph for the range from 0 cents to 2400 cents (two octaves).

```ts
import {type HarmonicEntropyOptions, precalculateRatios, harmonicEntropy} from 'harmonic-entropy';

const options: HarmonicEntropyOptions = {
  N: 10000,         // Maximum Benedetti height of rationals to consider (default)
  s: 0.01,          // Gaussian frequency deviation (default)
  a: 1,             // Rényi order (defaults to Shanon entropy)
  series: 'tenney', // Series of rationals to use (default)
  minCents: 0,      // Lower bound of tabulation (default)
  maxCents: 2400,   // Upper bound of tabulation (default)
  res: 1,           // Tabulation delta in cents (default)
  normalize: false, // Boolean flag to normalize the result by Hartley entropy (no normalization by default)
};

// Compute the set of rational numbers to consider.
const ratios = precalculateRatios(options);

// Compute the table of [cents, entropy] pairs. Entropy is measured in natural units.
const table = harmonicEntropy(options, ratios);

// This would be replaced by passing the table your favorite plotting library.
console.log(table);

/*
[
  [0, 2.465367706139234],
  [1, 2.4695232705775982],
  [2, 2.481975530994572],
  [3, 2.5026079139822173],
  [4, 2.5312627678840913],
  [5, 2.5676902925403278],
  ...
  [2399, 3.900196510219676],
  [2400, 3.898697709259859],
]
*/
```

Compute harmonic entropy of individual musical intervals.

```ts
import {EntropyCalculator} from 'harmonic-entropy';

const options = {maxCents: 1200};

// Preparing the internal tables takes some time.
const entropy = new EntropyCalculator(options);

// Evaluating entropy is fast.
const pureFifthEntropy = entropy.ofFraction(3/2);
console.log(pureFifthEntropy); // 4.121900707092091

// Tablulated values are linearly interpolated
const majorSixthEntropy = entropy.ofFraction('5/3');
console.log(majorSixthEntropy); // 4.42017406844399

// Intervals measured in cents are also supported in addition to frequency ratios.
const perfectFifthEntropy = entropy.ofCents(700);
console.log(perfectFifthEntropy); // 4.126342260377048

// Store internal tables for later.
const SERIALIZED = JSON.stringify(entropy);
```

Revive a serialized entropy calculator.

```ts
import {EntropyCalculator} from 'harmonic-entropy';

// Obtain SERIALIZED string from previous example.

// Revification is fast.
const entropy = JSON.parse(SERIALIZED, EntropyCalculator.reviver);

// Revived instance works like the original.
const perfectFifthEntropy = entropy.ofCents(700);
console.log(perfectFifthEntropy); // 4.126342260377048
```
