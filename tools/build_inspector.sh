echo Building JS Inspector
if [ ! -d "~inspect" ]; then
  echo Configuring JS Inspector
  mkdir ~inspect
  cd ~inspect && emconfigure ../../configure --disable-multithread --disable-runtime-cpu-detect --target=generic-gnu --enable-accounting --enable-inspection --enable-aom_highbitdepth --extra-cflags="-D_POSIX_SOURCE"
fi

cd ~inspect
emmake make -j 8
cp examples/inspect inspect.bc
emcc -O3 inspect.bc -o inspect.js -s TOTAL_MEMORY=134217728 -s MODULARIZE=1 -s EXPORT_NAME="'DecoderModule'" --post-js "../inspect-post.js" --memory-init-file 0
cp inspect.js ../inspect.js
