FROM schismdev/arch-linux-base:latest

RUN sudo -u aur yay --noconfirm -S wget

RUN mkdir -p /data/nwm_test_mesh && wget -qO- "https://www.dropbox.com/s/mjaxaqeggy721um/Gulf_Stream_develop.tar.gz?dl=1" | tar xz -C /data/nwm_test_mesh

ENV NWM_TEST_MESH /data/nwm_test_mesh/hgrid.gr3

RUN mkdir -p /root/.local/share/nwm \
	&& wget -P /root/.local/share/nwm -O NWM_channel_hydrofabric.tar.gz \
		https://www.nohrsc.noaa.gov/pub/staff/keicher/NWM_live/web/data_tools/NWM_channel_hydrofabric.tar.gz

RUN mkdir /opt/hostedtoolcache

ENV AGENT_TOOLSDIRECTORY /opt/hostedtoolcache
