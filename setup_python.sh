#!/bin/bash
# Quick start script for Python version of Q-POP Diffraction

echo "Q-POP Diffraction - Python Version Setup"
echo "========================================="
echo ""

# Check Python version
python_version=$(python3 --version 2>&1)
if [ $? -eq 0 ]; then
    echo "✓ Python found: $python_version"
else
    echo "✗ Python 3 not found. Please install Python 3.7 or higher."
    exit 1
fi

# Check if pip is available
if command -v pip3 &> /dev/null; then
    echo "✓ pip3 found"
    PIP_CMD="pip3"
elif command -v pip &> /dev/null; then
    echo "✓ pip found"
    PIP_CMD="pip"
else
    echo "✗ pip not found. Please install pip."
    exit 1
fi

# Install requirements
echo ""
echo "Installing Python dependencies..."
$PIP_CMD install -r requirements.txt

if [ $? -eq 0 ]; then
    echo "✓ Dependencies installed successfully"
else
    echo "✗ Failed to install dependencies"
    exit 1
fi

# Check if input files exist
echo ""
echo "Checking for required input files..."
required_files=("parameter.system.in" "parameter.atom.in")
optional_files=("phaseFra.in" "strucOrd.in" "displace.in" "region.in")

all_required_present=true
for file in "${required_files[@]}"; do
    if [ -f "$file" ]; then
        echo "✓ $file found"
    else
        echo "✗ $file NOT FOUND (required)"
        all_required_present=false
    fi
done

for file in "${optional_files[@]}"; do
    if [ -f "$file" ]; then
        echo "✓ $file found"
    else
        echo "  $file not found (optional - will use defaults)"
    fi
done

echo ""
if [ "$all_required_present" = true ]; then
    echo "Setup complete! You can now run:"
    echo "  python3 q-pop_diffraction.py"
else
    echo "Setup incomplete. Please provide the required input files."
    echo "You can find examples in the examples/ directory."
fi
