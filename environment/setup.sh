#!/bin/bash
# Setup ColabFold and PyTorch environments for the calibration kit
# Run: bash environment/setup.sh
# Then: source activate colabfold-env

set -e

echo "═══════════════════════════════════════════════════════════════"
echo "  Calibration Kit Environment Setup"
echo "═══════════════════════════════════════════════════════════════"
echo ""

# Check conda
if ! command -v conda &> /dev/null; then
    echo "✗ conda not found. Please install Miniconda or Anaconda."
    exit 1
fi

echo "✓ conda found: $(conda --version)"
echo ""

# Create ColabFold environment (JAX + CUDA)
echo "Setting up ColabFold environment (JAX, CUDA 12)..."
conda create -y -n colabfold-env \
    python=3.11 \
    pip \
    pytorch-cuda=12.1 \
    -c pytorch -c nvidia

echo "Activating colabfold-env..."
source activate colabfold-env

echo "Installing ColabFold and dependencies..."
pip install --upgrade pip
pip install colabfold[local]
pip install pandas scipy scikit-learn

echo ""
echo "Verifying environment..."
python -c "import colabfold; print(f'✓ colabfold {colabfold.__version__}')"
python -c "import torch; print(f'✓ torch {torch.__version__}')"
python -c "import pandas; print(f'✓ pandas {pandas.__version__}')"

echo ""
echo "═══════════════════════════════════════════════════════════════"
echo "✓ Setup complete!"
echo "═══════════════════════════════════════════════════════════════"
echo ""
echo "Next steps:"
echo "  1. Activate: source activate colabfold-env"
echo "  2. Run: python scripts/01_fold_benchmark.py --limit 5 (test with 5 pairs)"
echo "  3. See QUICKSTART.md for full workflow"
echo ""
