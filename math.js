const mapping = [
  'X_"D1"',
  'u_1',
  'X_"D2"',
  'u_2',
  'A_D',
  'K_"DD"',
  'K_"D1"',
  'K_"D2"',
  'K_"12"',
  'K_"1D"',
  'K_"11"',
  'K_"2D"',
  'K_"21"',
  'K_"22"',
];
const alphabet = [...'abcdefghijklmnopqrstuvwxyz'];

function encode(str) {
  let result = str;
  for (let i = 0; i < mapping.length; i++) {
    result = result.replaceAll(mapping[i], alphabet[i]);
  }
  return result;
}

function decode(str) {
  let result = str;
  for (let i = 0; i < mapping.length; i++) {
    result = result.replaceAll(alphabet[i], mapping[i]);
  }
  return result;
}

let radek1 = '(u_1*K_"11" + X_"D1"*u_1*K_"1D" + X_"D2"*u_2*K_"1D" + A_D*K_"1D" + u_2*K_"12")*(u_1)';
let result1 = 'a b^2 j + b^2 k + b c d j + i b d + b e j';

let radek2 = '(u_1*K_"D1"+X_"D1"*u_1*K_"DD"+X_"D2"*u_2*K_"DD"+A_D*K_"DD"+u_2*K_"D2")*(X_"D1"*u_1+X_"D2"*u_2+A_D)';
let result2 = 'e^2f+b^2a^2f+2ebaf+2bacdf+c^2d^2f+2ecdf+b^2ag+bcdg+ebg+badh+cd^2h+edh';

let radek3 = '(u_1*K_"21" + X_"D1"*u_1*K_"2D" + X_"D2"*u_2*K_"2D" + A_D*K_"2D" + u_2*K_"22")*(u_2)';
let result3 = 'mbd+lbad+lcd^2+eld+nd^2';


console.log(decode(result2));
