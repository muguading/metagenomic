CONDA_ROOT=$HOME/miniconda3
OUT=$HOME/pathogen-conda-packs

mkdir -p "$OUT"

$CONDA_ROOT/bin/conda install -n base -c conda-forge conda-pack -y
sudo apt-get update
sudo apt-get install -y zstd

pack_env() {
  src="$1"
  name="$2"
  tmp="$(mktemp -d)"
  echo "Packing $name"
  "$CONDA_ROOT/bin/conda-pack" -p "$src" -o "$tmp/$name.tar.gz" --force
  gzip -dc "$tmp/$name.tar.gz" | zstd -T0 -19 -o "$OUT/$name.tar.zst"
  rm -rf "$tmp"
}

for env_dir in "$CONDA_ROOT"/envs/*; do
  [ -d "$env_dir/conda-meta" ] || continue
  env_name="$(basename "$env_dir")"
  pack_env "$env_dir" "$env_name"
done
