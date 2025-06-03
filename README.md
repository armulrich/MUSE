# Fork of: MUSE 
Design and Analysis of Permanent Magnets for Stellarators

### Attempting to replicate finite mu section MUSE

My usage info: 
- using Python 3.12.9 
- Running on a HPC with hardware: ``srun -n 1 --pty --gres=gpu:rtx8000:1 -c 16 -t 2:00:00 --mem=100000 /bin/bash``

Instructions: 
1. If not already done create virtual environment and ``pip install --no-cache-dir -r requirements.txt``
2. Add the ```./Input/magtense_zot80_3d.csv``` to the project root
3. Open ``Attempt_to_Replicate_finite_mu_MUSE.ipynb``
4. Run the first cell which just runs ``!python mgrid_301.py --muea 1.05 --muoa 1.05``
   * Note: This is just one mu setting once this works correctly can cylce through many configurations of mu to try to replicate MUSE paper results. 
