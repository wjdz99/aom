#include <cstdio>

extern "C" int superfunc(void);

int main(int argc, const char* argv[]) {
  printf("Hello, super=%d!\n", superfunc());
  return 0;
}
