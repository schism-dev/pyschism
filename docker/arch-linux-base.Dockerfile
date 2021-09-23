from oblique/archlinux-yay

RUN pacman -Syu --noconfirm \
	&& sudo -u aur yay --noconfirm -S \
		netcdf-openmpi \
		pnetcdf-openmpi \
		netcdf-fortran-openmpi \
		udunits \
		python



RUN pacman -Qtdq | xargs -r pacman --noconfirm -Rcns
RUN rm -rf /home/aur/.cache