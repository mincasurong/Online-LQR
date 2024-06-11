# Online-LQR

This repository is dedicated to the implementation and study of Online Linear Quadratic Regulator (LQR), which combines traditional LQR with reinforcement learning techniques.

## Contents

- **MATLAB Files:**
  - **[CTLQR_PI_3x3_Sogang.m](CTLQR_PI_3x3_Sogang.m)**: Continuous-time LQR with Proportional-Integral control.
  - **[DTLQR_Q_3x3_Sogang.m](DTLQR_Q_3x3_Sogang.m)**: Discrete-time LQR with state-space representation.
  - **[DTLQR_Q_3x3_Sogang_Nonlin_ref.m](DTLQR_Q_3x3_Sogang_Nonlin_ref.m)**: Discrete-time LQR with nonlinear reference.
  - **[DTLQR_Q_3x3_Sogang_Nonlin_ref_PI.m](DTLQR_Q_3x3_Sogang_Nonlin_ref_PI.m)**: Discrete-time LQR with nonlinear reference and PI control.

- **Python Files:**
  - **[MPC_simple_pendulum.py](MPC_simple_pendulum.py)**: Model Predictive Control (MPC) implementation for a simple pendulum system.

## Getting Started

### Prerequisites

- MATLAB (for running the `.m` files)
- Python 3.6 or higher (for running the `.py` files)
- Required Python libraries: `numpy`, `scipy`, `matplotlib`

### Installation

1. Clone the repository:
    ```bash
    git clone https://github.com/mincasurong/Online-LQR.git
    cd Online-LQR
    ```

2. Install the required Python libraries:
    ```bash
    pip install numpy scipy matplotlib
    ```

## Usage

### Running MATLAB Files

Open the `.m` files in MATLAB and run them to see the results of various LQR implementations.

### Running Python Files

To run the Model Predictive Control example for a simple pendulum:
```bash
python MPC_simple_pendulum.py
