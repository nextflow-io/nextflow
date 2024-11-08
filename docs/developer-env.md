(developer-env)=

# Environment setup

Setting up a Nextflow development environment is a prerequisite for creating, testing, and optimizing data analysis pipelines. The steps below outline recommended tools for setting up an optimal Nextflow development environment.

:::{note}
If you are using a Windows computer, you first need to install and configure the Windows Subsystem for Linux (WSL). See {ref}`wsl` for installation instructions.
:::

(vscode-install)=

## VS Code

Installing an Integrated Development Environment (IDE) is an essential step for setting up your environment and provide a user-friendly interface for writing, editing, and managing code.

Visual Studio Code (VS Code) is a popular lightweight IDE that is known for its versatility and extensibility. It offers features like syntax highlighting, intelligent code completion, and integrated debugging tools for various programming languages. VS Code supports Windows, macOS, and Linux, and is a good choice for both new and experienced Nextflow developers.

````{tabs}

```{group-tab} Windows

To install VS Code on Windows:

1. Visit the [VS Code](https://code.visualstudio.com/download) website.
1. Download VS Code installer for Windows.
1. Double-click the installer executable (`.exe`) file and follow the step-by-step setup installation wizard.

```

```{group-tab} macOS

To install VS Code on macOS:

1. Visit the [VS Code](https://code.visualstudio.com/download) website.
1. Download VS Code for macOS.
1. Drag the `Visual Studio Code.app` application to the Applications folder to make it available in the macOS Launchpad.

```

```{group-tab} Linux

To install VS Code on Linux Debian/Ubuntu distributions:

1. Visit the [VS Code](https://code.visualstudio.com/download) website.
1. Download VS Code Linux Debian/Ubuntu (`.deb`) distribution.
1. Open a new terminal window.
1. Navigate to the folder containing your VS Code download.
1. Run `sudo apt install ./<file>.deb`, replacing `<file>` with the full file name.
   
   :::{note}
   If you're using an older Linux distribution, run `sudo dpkg -i <file>.deb` to install VS Code and `sudo apt-get install -f` to install dependencies.
   :::
   
See [Linux installation](https://code.visualstudio.com/docs/setup/linux#_installation) for information about installing VS Code on other distributions.

```

````

## Extensions

Extensions are a key feature of IDE's and allow you to customize your development environment by adding support for various programming languages, tools, and features. The [VS Code Marketplace](https://marketplace.visualstudio.com/vscode) offers thousands of extensions that can enhance your productivity and tailor the editor to your specific needs.


### Nextflow

The VS Code [Nextflow extension](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow) adds Nextflow language support to the editor. The Nextflow extension enhances development with:

- Diagnostics
- Hover hints
- Code navigation
- Code completion
- Formatting
- Renaming
- Parameter schemas
- DAG previews

See {ref}`vs-code-page` for more information about the Nextflow extension.

````{tabs}

```{group-tab} Windows

To install the Nextflow VS Code extension on Windows:

1. Open VS Code.
1. Open the VS Code Extensions view in the left-hand menu.
1. Search for Nextflow.
1. Select **Install**.

```

```{group-tab} macOS

To install the Nextflow VS Code extension on macOS:

1. Open VS Code.
1. Open the VS Code Extensions view in the left-hand menu.
1. Search for Nextflow.
1. Select **Install**.

```

```{group-tab} Linux

To install the Nextflow VS Code extension on Linux Debian/Ubuntu distributions:

1. Open VS Code.
1. Open the VS Code Extensions view in the left-hand menu.
1. Search for Nextflow.
1. Select **Install**.

```

````

(remote-development-ext)=

### Remote Development

The [Remote Development extension pack](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.vscode-remote-extensionpack) enables you to run WSL, SSH, or a development container for editing and debugging with the full set of VS Code features.

The Remote Development extension pack includes four extensions:

[Remote - SSH](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-ssh) - Work with source code in any location by opening folders on a remote machine/VM using SSH.
[Remote - Tunnels](https://marketplace.visualstudio.com/items?itemName=ms-vscode.remote-server) - Work with source code in any location by opening folders on a remote machine/VM using a VS Code Tunnel (rather than SSH).
[Dev Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) - Work with a separate toolchain or container based application by opening any folder mounted into or inside a container.
[WSL](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-wsl) - Get a Linux-powered development experience from the comfort of Windows by opening any folder in the Windows Subsystem for Linux.

:::{note}
The VS Code Remote Development extension pack is required if you are developing using remote servers, Windows Subsystem for Linux, or Development Containers.
:::


````{tabs}

```{group-tab} Windows

To install the Remote Development extension pack on Windows:

1. Open VS Code.
1. Open the VS Code Extensions view in the left-hand menu.
1. Search for Nextflow.
1. Select **Install**.

```

```{group-tab} macOS

To install the Remote Development extension pack on macOS:

1. Open VS Code.
1. Open the VS Code Extensions view in the left-hand menu.
1. Search for Nextflow.
1. Select **Install**.

```

```{group-tab} Linux

To install the Remote Development extension pack on Linux Debian/Ubuntu distributions:

1. Open VS Code.
1. Open the VS Code Extensions view in the left-hand menu.
1. Search for Nextflow.
1. Select **Install**.

```

````

(docker-desktop)=

## Docker

Containerization enables the creation of self-contained and fully reproducible computational pipelines by bundling a script's binary dependencies into a standardized and portable format. Containers can be executed on any platform that supports a container runtime and ensures consistency across different environments.

Docker is an open-source platform that simplifies application development, deployment, and execution by packaging applications and their dependencies into containers. Docker Desktop provides a Graphical User Interface (GUI) for managing Docker containers. Installing Docker Desktop is a straightforward process that allows you to create, deploy, and manage applications within containers.

<!---
Configure your environment to support the container technologies you want to use. 
--->

````{tabs}

```{group-tab} Windows

To install Docker Desktop on Windows:

1. Visit the [Install Docker Desktop on Windows](https://docs.docker.com/desktop/setup/install/windows-install/) page.
1. Download the installer using the download button at the top of the page, or from the release notes.
1. Double-click Docker Desktop `Installer.exe` to run the installer. By default, Docker Desktop is installed at `C:\Program Files\Docker\Docker`.
1. When prompted, ensure the Use WSL 2 instead of Hyper-V option on the Configuration page is selected, or not, depending on your choice of backend.

    :::{note}
    You won't be able to select which backend to use if your system only supports one of the two options.
    :::

1. Follow the instructions on the installation wizard to authorize the installer and proceed with the install.
1. When the installation is complete, select **Close**.
1. Start Docker Desktop.
1. Review the Docker Subscription Service Agreement and, if you agree, select **Accept** to continue.

    :::{note}
    Docker Desktop won't run if you do not agree to the terms. You can choose to accept the terms at a later date by opening Docker Desktop.
    :::

1. Docker Desktop starts after you accept the terms.

```

```{group-tab} macOS

To install Docker Desktop on macOS:

1. Visit the [Install Docker Desktop on Mac](https://docs.docker.com/desktop/install/mac-install/) page.
1. Download the installer for your chip type using the download buttons at the top of the page.
1. Double-click `Docker.dmg` to open the installer.
1. Drag the Docker icon to the **Applications** folder to make it available in the macOS Launchpad.
1. Double-click **Docker.app** in the **Applications** folder to start Docker.
1. Review the Docker Subscription Service Agreement and, if you agree, select **Accept** to continue.
1. From the installation window, select **Use recommended settings (Requires password)**.

    :::{note} Recommended settings lets Docker Desktop automatically set the necessary configuration settings. Advanced settings allow you to set the location of the Docker CLI tools either in the system or user directory, enable the default Docker socket, and enable privileged port mapping. See [Settings](https://docs.docker.com/desktop/settings/#advanced), for more information and how to set the location of the Docker CLI tools.
    :::

1. Select **Finish**. If you have applied any of the previous configurations that require a password, enter your password to confirm your choice.

```

```{group-tab} Linux

To install Docker Desktop on Linux Debian/Ubuntu distributions:

1. Set up Docker's package repository. See step one of [Install using the `apt` repository](https://docs.docker.com/engine/install/ubuntu/#install-using-the-repository).
1. Download the latest Debian/Ubuntu (`.deb`) distribution.
1. In your terminal, run `sudo apt-get install ./docker-desktop-amd64.deb`

    :::{note}
    By default, Docker Desktop is installed at `/opt/docker-desktop`.
    :::

1. Double-click **Docker Desktop** in your Applications menu to start Docker.
1. Review the Docker Subscription Service Agreement and, if you agree, select **Accept** to continue.
1. From the installation window, select **Use recommended settings (Requires password)**. Docker Desktop starts after you accept the terms.

    :::{note}
    Docker Desktop won't run if you do not agree to the terms. You can choose to accept the terms at a later date by opening Docker Desktop.
    :::

```

````

Nextflow supports multiple container technologies (e.g., Singularity and Podman) allowing you to choose the one that best fits your needs. See {ref}`container-page` for more information about other supported container engines.

<!---
## Conda

Conda is an open-source package and environment manager that simplifies installing and configuring complex software across platforms. Nextflow supports Conda, enabling the use of Conda recipes and environment files to configure workflow dependencies.

:::{note}
Conda environments can lead to inconsistencies across systems due to dependency resolution and OS variations.
:::

The preferred method for installing Conda is through Miniconda, a lightweight version of Anaconda that includes Conda and its dependencies.

````{tabs}

```{group-tab} Windows

```

```{group-tab} macOS

To install Conda on macOS:

1. Visit the [Miniconda](https://docs.anaconda.com/miniconda/#miniconda) website.
1. Download the latest version of the `.pkg` Miniconda installer. 
1. Double-click the `.pkg` file.
1. Follow the step-by-step setup installation instructions.
1. When the installation finishes, open a new terminal window and run `conda list` to verify Conda was installed correctly.

See [Quick command line install](https://docs.anaconda.com/miniconda/#quick-command-line-install) for command line installation instructions.

```

```{group-tab} Linux

To install Conda on Linux Debian/Ubuntu distributions:

1. Visit the [Miniconda](https://docs.anaconda.com/miniconda/#miniconda) website.
1. Download the latest version of the `.sh` Miniconda installer.
1. In your terminal, run `bash <file name>.sh`, replacing `<file name>` with the installer file name.
1. Follow the step-by-step setup installation prompts.
1. When the installation finishes, open a new terminal window and run `conda list` to verify Conda was installed correctly.

```

````
--->

## Git

Git provides powerful version control that helps track code changes. Git operates locally, meaning you don't need an internet connection to track changes, but it can also be used with remote platforms like GitHub, GitLab, or Bitbucket for collaborative development.

Nextflow seamlessly integrates with Git for source code management providers for managing pipelines as version-controlled Git repositories.

````{tabs}

```{group-tab} Windows

Git is already installed on most WSL distributions. You can check if it is already installed by running `git version`.

To install the latest stable Git version on Linux Debian/Ubuntu distributions:

1. Open a terminal window and run `sudo apt-get install git-all`.
1. Once complete, run `git version` to verify Git was installed.

See [git-scm documentation](https://git-scm.com/downloads/linux) for more information about installing Git on other Linux distributions.

```

```{group-tab} macOS

Git installed is already installed on new versions of macOS. You can activate it through the terminal running `git version`. If Git is not installed, you can install the latest version of Git using several methods:

To install Git on macOS with [Homebrew](https://docs.brew.sh/):

1. Open a terminal window and run `brew install git`.

   :::{note}
   You must have Homebrew installed. See [Homebrew installation](https://docs.brew.sh/Installation) for instructions.
   :::

1. Once complete, run `git version` to verify Git was installed.

To install Git on macOS with [Xcode](https://developer.apple.com/xcode/):

1. Open the App Store on your Mac.
1. Sign in to your Apple Account.
1. Search for Xcode.
1. Select **Install**.
1. Once complete, open a new terminal window and run `git version` to verify Git was installed.

```

```{group-tab} Linux

Git is already installed on most on most Linux Debian/Ubuntu distributions.

To install the latest stable Git version on Linux Debian/Ubuntu distributions:

1. Open a terminal window and run `sudo apt-get install git-all`.
1. Once complete, run `git version` to verify Git was installed.

See [git-scm documentation](https://git-scm.com/downloads/linux) for more information about installing Git on other Linux distributions.

```

````

(wsl)=

## Windows Subsystem for Linux

Developers can access the power of both Windows and Linux on a Windows machine. The Windows Subsystem for Linux (WSL) lets developers install a Linux distribution and use Linux applications, utilities, and Bash command-line tools directly on Windows without the overhead of a virtual machine or dual-boot setup.

WSL is an optional feature on Windows 10 version 2004 and higher (Build 19041 and higher) or Windows 11. You can enable it through PowerShell or Windows Command Prompt. The steps below outline the recommended setup.

<!---
### Windows Features dialog

To enable WSL on Windows using Windows Features dialog:

1. In the Windows search bar, enter 'features' to bring up the **Turn Windows Features on and off** dialog.
1. Scroll down and check **Windows Subsystem for Linux**.
1. Select **OK** and restart Windows.
1. After restarting Windows, check that you have WSL enabled by opening a Command Prompt or PowerShell and typing `wsl`.

To install Ubuntu on WSL:

1. Go to **Start Button > Microsoft Store**.
1. Enter 'Linux' into the search field, then click **Run Linux on Windows**.
1. Select the latest Ubuntu distribution.
1. Select **Get** and wait for Windows to download and install Ubuntu.
1. When it’s finished, select **Launch**.
1. A terminal window will appear. Wait for Ubuntu to finish installing, then create a new Linux username and password when prompted.

### PowerShell or Windows Command Prompt
--->

To enable WSL on Windows using Powershell or Windows Command Prompt:

1. Open PowerShell or Windows Command Prompt in administrator mode by right-clicking and selecting **Run as administrator**.
1. Run `wsl --install`.

    :::{note}
    This command will enable the features necessary to run WSL and install the Ubuntu distribution.
    :::

1. When prompted, restart Windows.
1. After restarting Windows, open the Ubuntu distribution using the **Start** menu and create a new Linux **User Name** and **Password** when prompted.

    :::{note}
    The **User Name** and **Password** is specific to each Linux distribution that you install and has no bearing on your Windows user name.
    :::

See [Set up a WSL development environment](https://learn.microsoft.com/en-us/windows/wsl/setup/environment) for more information.

## Development Containers

Development Containers (Dev Containers), are Docker containers that are specifically configured to provide a fully featured development environment. They can be used to run an application, to separate tools, libraries, or runtimes needed for working with a codebase. Dev Containers can be run locally or remotely, in a private or public cloud, and in a variety of supporting tools and editors.

VS code and Docker are required to create and manage your Dev Containers. See {ref}`vscode-install` and {ref}`docker-desktop` for installation instructions.

### Development Containers extension

The VS Code Dev Containers extension lets you use a container as a full-featured development environment. It allows you to open any folder inside (or mounted into) a container and take advantage of VS Code's full feature set.

The Dev Containers extension is included as a part of the [Remote Development extension pack](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.vscode-remote-extensionpack). See {ref}`remote-development-ext` for installation instructions.

### Create and run a dev container

A `devcontainer.json` file in your project directory instructs VS Code how to access, create, and configure a development container. It can be used to run an application or to separate tools, libraries, or runtimes.

The `devcontainer.json` can be used to:

- Install additional tools in the container.
- Automatically install extensions.
- Forward or publish additional ports.
- Set runtime arguments.
- Reuse or extend your existing Docker Compose setup.
- Add more Advanced container configuration.

To create a dev container with an existing image:

1. Create a Dev Container `.json` under `.devcontainer/devcontainer.json` in the root of your project.
1. Add an image with the Nextflow tooling, and VS Code customizations, to the `.json` file. For example:

    ```json
    {
        "name": "Nextflow Dev Container",
        "image": "nfcore/gitpod:latest",
        "remoteUser": "vscode",
        "runArgs": ["--privileged"],

        // Configure tool-specific properties.
        "customizations": {
            // Configure properties specific to VS Code.
            "vscode": {
                // Set *default* container specific settings.json values on container create.
                "settings": {
                    "python.defaultInterpreterPath": "/opt/conda/bin/python"
                },

                // Add the IDs of extensions you want installed when the container is created.
                "extensions": ["ms-python.python", "ms-python.vscode-pylance", "nf-core.nf-core-extensionpack", "nextflow.nextflow"]
            }
        }
    }
    ```

    :::{note}
    Instead of using a prebuilt image, a custom Dockerfile may also live in the `.devcontainer` folder. You can replace the image property in `devcontainer.json` with dockerfile and utilize the custom container. See [Create a Dev Container](https://code.visualstudio.com/docs/devcontainers/create-dev-container) for more information.
    :::

1. Enter **Dev Containers: Reopen in Container** in the VS Code Command Palette and reopen your project. You should now see the name of the container ("Nextflow Dev Container" if using the above example) in the bottom left corner of VS Code.

:::{note}
Dev Containers can also be used by GitHub Codespaces in VS Code or the browser. See [GitHub Codespaces](https://code.visualstudio.com/docs/remote/codespaces) for more information.
:::
