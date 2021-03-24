import matplotlib.pyplot as plt
import numpy as np
import os
import imageio
import shutil



def read(name, ext='', directory=None, log_dens=False):
	if directory is None:
		directory = rf'{os.getcwd()}/{name}{ext}'

	print(directory)

	try: os.mkdir(directory + "/frames")
	except: pass
	try: os.mkdir(directory + f"/frames/{name}{ext}")
	except: pass
	matrices = []
	eigenval = []
	energies = []
	true_energies = []
	iterations = 0
	with open(directory + f"/hf_{name}.out", 'r') as f:
		lines = f.readlines()
		recorddens = False
		recorden = False
		for line in lines:
			if line.strip().startswith('Total energy'):
				parts = line.split(':')
				true_energies.append(float(parts[-1].strip()))

			if line.strip().startswith('Energy:'):
				parts = line.split(':')
				energies.append(float(parts[-1].strip()))

			if 'end density' in line and log_dens:
				recorddens = False
				continue
			elif 'begin density' in line and log_dens:
				recorddens = True
				matrices.append([])
				counter = 0
				continue

			if 'end eigenvalues' in line:
				recorden = False
				continue
			elif 'begin eigenvalues' in line:
				recorden = True
				eigenval.append([])
				counter = 0
				iterations += 1
				continue

			if recorddens and log_dens:
				parts = line.split()
				matrices[-1].append([])
				for k in range(len(parts)):
					matrices[-1][counter].append(float(parts[k]))
				counter += 1

			if recorden:
				eigenval[-1].append(float(line))
				counter += 1

	if log_dens:
		frames = []
		print(f"Found density matrices for {len(matrices)} SCF iterations...")
		for i in range(len(matrices)//2):
			print(f'{i:>4}: HF energy = {energies[i]:7.9f} hartree')
			fig = plt.figure()
			fig.suptitle(f"Iteration {i}, Energy = {energies[i]:7.9f}")
			ax1 = fig.add_subplot(121)
			ax1.imshow(matrices[2*i], cmap='viridis')
			ax1.title.set_text("Real part of density")
			ax2 = fig.add_subplot(122)
			ax2.imshow(matrices[2*i+1], cmap='viridis')
			ax2.title.set_text("Imag part of density")
			heatmap = plt.pcolor(matrices[2*i])
			plt.colorbar(heatmap)
			fig.savefig(directory + "/frames" + f'/{i}.png')
			frames.append(directory + "/frames" + f'/{i}.png')

			plt.close()

		images = [imageio.imread(f) for f in frames]
		imageio.mimsave(directory + f'/anim.gif', images, fps=5)

	fig = plt.figure()
	plt.plot(range(iterations), energies, label='HF Energy')
	plt.plot(range(iterations), [e[0] for e in eigenval], label='Lowest eigenval')
	plt.legend()
	plt.title("Per Iteration")
	plt.xlabel('Iteration')
	plt.ylabel('Energy (hartree)')
	fig.savefig(directory + f'/energy.jpg')

	plt.close()

	return energies, true_energies


if __name__ == '__main__':
	name = 'h2'
	ext = '_UHF'

	# read(name, ext, directory=rf"{os.getcwd()}/{name}{ext}")


	try: os.mkdir(f"{os.getcwd()}/h2_r")
	except: pass

	energies = []
	true_energies = []
	r_range = np.linspace(0.4, 3, 15)
	for i, r in enumerate(r_range):
		try: os.mkdir(f"{os.getcwd()}/h2_r/{i}")
		except: pass

		shutil.copy(f"{os.getcwd()}/h2_r/hf.inp", f"{os.getcwd()}/h2_r/{i}/hf.inp")
		shutil.copy(f"{os.getcwd()}/h2_r/settings", f"{os.getcwd()}/h2_r/{i}/settings")

		with open(f"{os.getcwd()}/h2_r/{i}/h2.xyz", 'w+') as f:
			f.write(f"2 \nhydrogen with bond stretched to twice the normal length \nH     0.0 0.0 0.0 \nH     0.0 0.0 {r} \n")

		os.system(f'pam --inp={os.getcwd()}/h2_r/{i}/hf.inp --mol={os.getcwd()}/h2_r/{i}/h2.xyz --copy={os.getcwd()}/h2_r/{i}/settings')
		#move files
		shutil.move(f"{os.getcwd()}/hf_h2.out", f"{os.getcwd()}/h2_r/{i}/hf_h2.out")
		shutil.move(f"{os.getcwd()}/hf_h2.tgz", f"{os.getcwd()}/h2_r/{i}/hf_h2.tgz")

		res = read('h2', f'_{i}', f'{os.getcwd()}/h2_r/{i}')

		energies.append(res[0][-1])
		true_energies.append(res[1][-1])


	fig = plt.figure()
	plt.plot(r_range, energies, label='Yuman')
	plt.plot(r_range, true_energies, label='Dirac')
	plt.legend()
	plt.title("PES bond distance")
	plt.xlabel('r (Angstrom)')
	plt.ylabel('HF Energy (hartree)')
	fig.savefig(f"{os.getcwd()}/PES.jpg")

	for i, e, te in zip(range(len(energies)), energies, true_energies):
		print(f"Step {i:>3}: E = {e:.6f} E_dirac = {te:.6f}")

	print('\a')