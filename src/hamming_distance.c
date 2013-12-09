///////////////////////////////////////////////////////////////////////
// Copied from phash (http://www.phash.org/): src/pHash.cpp to avoid installing all the library and its dependences.
// We would like to thank Evan Klinger & David Starkweather for their library
typedef unsigned long long ulong64;
int ph_hamming_distance_64bit(const ulong64 hash1,const ulong64 hash2){
  ulong64 x = hash1^hash2;
  const ulong64 m1  = 0x5555555555555555ULL;
  const ulong64 m2  = 0x3333333333333333ULL;
  const ulong64 h01 = 0x0101010101010101ULL;
  const ulong64 m4  = 0x0f0f0f0f0f0f0f0fULL;
  x -= (x >> 1) & m1;
  x = (x & m2) + ((x >> 2) & m2);
  x = (x + (x >> 4)) & m4;
  return (x * h01)>>56;
}
///////////////////////////////////////////////////////////////////////

ulong64 hatol(char *s) {
  if (s == 0) return 0;
  if (*s && *s == '0' && *(s + 1) && (*(s+1) == 'x' || *(s+1) == 'X'))
    s += 2;

  ulong64 d,l;
  for(l = 0; *s ; ++s){
    if (*s >= '0' && *s <= '9')
      d = *s - '0';
    else if (*s >= 'a' && *s <= 'f')
      d = *s - 'a' + 10;
    else if (*s >= 'A' && *s <= 'F')
      d = *s - 'A' + 10;
    else
      break;
    l = (l << 4) + d;
  }

  return (l);
}
int ph_hamming_distance(char* str_hash1, char* str_hash2){
  const ulong64 hash1 = hatol(str_hash1);
  const ulong64 hash2 = hatol(str_hash2);
  return ph_hamming_distance_64bit(hash1,hash2);
}
