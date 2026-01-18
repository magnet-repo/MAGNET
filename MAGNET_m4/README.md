# MAGNET_m4

This folder contains the **MAGNET_m4** implementation targeting Nordic Semiconductor platforms based on an **ARM Cortex-M4** CPU. The implementation is built and flashed using **Zephyr** via the **nRF Connect SDK (NCS) v3.1.1**.

The implementation was evaluated on the **nRF52840 Development Kit**. Other Nordic nRF development kits equipped with an ARM Cortex-M4 CPU should also be compatible.

---

## Requirements

- Linux (recommended)
- Nordic **nrfutil**
- nRF Connect SDK (**NCS**) v3.1.1
- Supported Nordic development kit (for example, `nrf52840dk`)

---

## 1. Installing `nrfutil`

The `nrfutil` command-line tool is required to install and manage the nRF Connect SDK toolchain.

### 1.1 Download `nrfutil`

Download the latest standalone `nrfutil` binary from Nordic Semiconductor’s official GitHub releases page:

[https://github.com/NordicSemiconductor/nrfutil/releases]

### 1.2 Extract and install

Extract the archive and move the `nrfutil` binary to a directory included in your `PATH`:

```
tar -xvf nrfutil-linux-x86_64.tar.gz
chmod +x nrfutil
sudo mv nrfutil /usr/local/bin/
```
### 1.3 Verify the installation

```
nrfutil --version
```

### 1.4 Launching the nRF Connect SDK Toolchain

Once nrfutil is installed, launch a shell with the toolchain corresponding to nRF Connect SDK v3.1.1:

```
nrfutil sdk-manager toolchain launch --ncs-version v3.1.1 --shell
```

### 1.5 Fetching the nRF Connect SDK Workspace

Initialize a new workspace using the official NCS manifest:

```
west init -m https://github.com/nrfconnect/sdk-nrf --mr v3.1.1 v3.1.1
cd v3.1.1
west update
```

This step downloads Zephyr and all required NCS components.

## 2. Building MAGNET_m4

From the Zephyr directory within the NCS workspace, build the application for the nRF52840 development kit:

```
cd ncs/v3.1.1/zephyr
west build -b nrf52840dk/nrf52840 MAGNET_m4
```

### 2.1 Flashing the Firmware

Connect the development kit to the host machine via USB and flash the firmware:

```
west flash
```

### 2.2 Viewing Output (Optional)

Program output can be monitored using the Serial Terminal provided by the nRF Connect for Desktop application or any compatible serial terminal. Ensure that the correct serial device is selected.




























