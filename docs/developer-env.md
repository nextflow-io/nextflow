(devenv-page)=

# Environment setup

Setting up a Nextflow development environment is a prerequisite for creating, testing, and optimizing data analysis pipelines. The steps below outline recommended tools for setting up an optimal Nextflow development environment.

Nextflow must be installed separately. See {ref}`install-page` for Nextflow installation instructions.

:::{note}
If you are using a Windows computer, first install and configure the Windows Subsystem for Linux (WSL). See {ref}`devenv-wsl` for installation instructions.
:::

(devenv-vscode)=

## VS Code

An Integrated Development Environment (IDE) provides a user-friendly interface for writing, editing, and managing code. Installing one is an essential step for setting up your environment. 

Visual Studio Code (VS Code) is a popular lightweight IDE known for its versatility and extensibility. It offers features like syntax highlighting, intelligent code completion, and integrated debugging tools for various programming languages. VS Code supports Windows, macOS, and Linux, and is a good choice for both new and experienced Nextflow developers.

````{tabs}

```{group-tab} Windows

To install VS Code on Windows:

1. Visit the [VS Code](https://code.visualstudio.com/download) website.
1. Download VS Code for Windows.
1. Double-click the installer executable (`.exe`) file and follow the set up steps.

```

```{group-tab} macOS

To install VS Code on macOS:

1. Visit the [VS Code](https://code.visualstudio.com/download) website.
1. Download VS Code for macOS.
1. Drag the `Visual Studio Code.app` application to the Applications folder to add it to the macOS Launchpad.

```

```{group-tab} Linux

To install VS Code on Linux Debian/Ubuntu distributions:

1. Visit the [VS Code](https://code.visualstudio.com/download) website.
1. Download the VS Code Linux Debian/Ubuntu (`.deb`) distribution.
1. Open a new terminal window.
1. Navigate to the folder where you downloaded VS Code.
1. Run `sudo apt install ./<file>.deb`, replacing `<file>` with the full file name.
   
   :::{note}
   If you're using an older Linux distribution, run `sudo dpkg -i <file>.deb` to install VS Code and `sudo apt-get install -f` to install dependencies.
   :::
   
See [Linux installation](https://code.visualstudio.com/docs/setup/linux#_installation) for information about installing VS Code on other distributions.

```

````

## Extensions

Extensions are a key feature of IDEs and allow you to customize your development environment by adding support for various programming languages, tools, and features. The [VS Code Marketplace](https://marketplace.visualstudio.com/vscode) offers thousands of extensions that can enhance your productivity and tailor the editor to your specific needs. Popular VS Code extensions for Nextflow developers are listed below:

**Nextflow**

The VS Code [Nextflow extension](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow) adds Nextflow language support to the editor. The Nextflow extension enhances development with:

- Diagnostics
- Hover hints
- Code navigation
- Code completion
- Formatting
- Renaming
- Parameter schemas
- DAG previews

See {ref}`vscode-page` for more information about the Nextflow extension features and how it enforces the Nextflow syntax.

**nf-core**

The [nf-core extension pack](https://marketplace.visualstudio.com/items?itemName=nf-core.nf-core-extensionpack) adds a selection of tools that help develop with nf-core, a community effort to collect a curated set of analysis pipelines built using Nextflow.

The nf-core extension pack includes several useful extensions. For example, [Code Spell Checker](https://marketplace.visualstudio.com/items?itemName=streetsidesoftware.code-spell-checker), [Prettier](https://marketplace.visualstudio.com/items?itemName=esbenp.prettier-vscode), [Todo Tree](https://marketplace.visualstudio.com/items?itemName=Gruntfuggly.todo-tree), and [Markdown Extended](https://marketplace.visualstudio.com/items?itemName=jebbs.markdown-extended). See [nf-core extension pack](https://marketplace.visualstudio.com/items?itemName=nf-core.nf-core-extensionpack) for more information about the tools included in the nf-core extension pack.

(devenv-remote)=

**Remote development**

The [Remote Development extension pack](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.vscode-remote-extensionpack) enables you to run WSL, SSH, or a development container for editing and debugging with the full set of VS Code features.

The pack includes the [Remote - SSH](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-ssh), [Remote - Tunnels](https://marketplace.visualstudio.com/items?itemName=ms-vscode.remote-server), [Dev Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers), and [WSL](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-wsl) extensions. See [Remote Development](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.vscode-remote-extensionpack) for more information about the tools included in the remote development extension pack.

:::{note}
The Remote Development extension pack is required if you are developing using remote servers, Windows Subsystem for Linux, or Development Containers.
:::

Installing VS Code extensions requires just a few clicks in the Extensions Marketplace. 

````{tabs}

```{group-tab} Windows

To install a VS Code extension on Windows:

1. Open VS Code.
1. Open the **Extensions** view in the left-hand menu.
1. Search for the extension.
1. Select **Install**.

```

```{group-tab} macOS

To install a VS Code extension on macOS:

1. Open VS Code.
1. Open the **Extensions** view in the left-hand menu.
1. Search for the extension.
1. Select **Install**.

```

```{group-tab} Linux

To install a VS Code extension on Linux Debian/Ubuntu distributions:

1. Open VS Code.
1. Open the **Extensions** view in the left-hand menu.
1. Search for the extension.
1. Select **Install**.

```

````

(devenv-docker)=

## Docker

Docker is an open-source platform that simplifies application development, deployment, and execution by packaging applications and their dependencies into containers. Containerization enables the creation of self-contained and fully reproducible computational pipelines by bundling a script's binary dependencies into a standardized and portable format. Containers can be executed on any platform that supports a container runtime and ensures consistency across different environments.

Docker Desktop provides a Graphical User Interface (GUI) for managing Docker containers. Installing Docker Desktop is a straightforward process that allows you to create, deploy, and manage applications within containers.

````{tabs}

```{group-tab} Windows

To install Docker Desktop on Windows:

1. Go to [Install Docker Desktop on Windows](https://docs.docker.com/desktop/setup/install/windows-install/).
1. Download the installer.
1. Double-click Docker Desktop `Installer.exe` to run the installer. By default, Docker Desktop is installed at `C:\Program Files\Docker\Docker`.
1. Depending on your choice of backend, select the **Use WSL 2 instead of Hyper-V** option on the Configuration page.

    :::{note}
    You won't be able to select which backend to use if your system only supports one of the two options.
    :::

1. Follow the instructions on the installation wizard to authorize the installer and proceed with the install.
1. When the installation is complete, select **Close**.
1. Start Docker Desktop.
1. Review the Docker Subscription Service Agreement and select **Accept** to continue.

    :::{note}
    Docker Desktop won't run if you do not agree to the terms. You can choose to accept the terms at a later date by opening Docker Desktop.
    :::

1. Docker Desktop starts after you accept the terms.

```

```{group-tab} macOS

To install Docker Desktop on macOS:

1. Go to [Install Docker Desktop on Mac](https://docs.docker.com/desktop/install/mac-install/).
1. Download the installer for your chip type.
1. Double-click `Docker.dmg` to open the installer.
1. Drag the Docker icon to the **Applications** folder to add it to the macOS Launchpad.
1. Double-click **Docker.app** in the **Applications** folder to start Docker.
1. Review the Docker Subscription Service Agreement and select **Accept** to continue.

    :::{note}
    Docker Desktop won't run if you do not agree to the terms. You can choose to accept the terms at a later date by opening Docker Desktop.
    :::
    
1. From the installation window, select **Use recommended settings (Requires password)**.

    :::{note} The **recommended settings** allow Docker Desktop to automatically set the necessary configuration settings. Advanced settings allow you to set the location of the Docker CLI tools either in the system or user directory, enable the default Docker socket, and enable privileged port mapping. See [Settings](https://docs.docker.com/desktop/settings/#advanced), for more information and how to set the location of the Docker CLI tools.
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
1. Review the Docker Subscription Service Agreement and select **Accept** to continue.
1. From the installation window, select **Use recommended settings (Requires password)**. Docker Desktop starts after you accept the terms.

    :::{note}
    Docker Desktop won't run if you do not agree to the terms. You can choose to accept the terms at a later date by opening Docker Desktop.
    :::

```

````

Nextflow supports multiple container technologies (e.g., Singularity and Podman) so you can choose the one that best fits your needs. See {ref}`container-page` for more information about other supported container engines.

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

To install Git on macOS with [Homebrew](https://docs.brew.sh/):

1. Open a terminal window and run `brew install git`.

   :::{note}
   You must have Homebrew installed. See [Homebrew installation](https://docs.brew.sh/Installation) for instructions.
   :::

1. Once complete, run `git version` to verify Git was installed.

To install Git on macOS with [Xcode](https://developer.apple.com/xcode/):

1. Open the App Store on your Mac.
1. Search for Xcode.
1. Select **Install**.
1. Once complete, open a new terminal window and run `git version` to verify Git was installed.

```

```{group-tab} Linux

Git is already installed on most Linux Debian/Ubuntu distributions.

To install the latest stable Git version on Linux Debian/Ubuntu distributions:

1. Open a terminal window and run `sudo apt-get install git-all`.
1. Once complete, run `git version` to verify Git was installed.

See [git-scm documentation](https://git-scm.com/downloads/linux) for more information about installing Git on other Linux distributions.

```

````

(devenv-wsl)=

## Windows Subsystem for Linux

Developers can access the power of both Windows and Linux on a Windows machine. The Windows Subsystem for Linux (WSL) lets developers install a Linux distribution and use Linux applications, utilities, and Bash command-line tools directly on Windows without the overhead of a virtual machine or dual-boot setup.

WSL is an optional feature on Windows 10 version 2004 and higher (Build 19041 and higher) or Windows 11. You can enable it through PowerShell or Windows Command Prompt. The steps below outline the recommended setup.

To enable WSL on Windows using Powershell or Windows Command Prompt:

1. Right-click and select **Run as administrator** to use PowerShell or Windows Command Prompt in administrator mode.
1. Run `wsl --install`.

    :::{note}
    This command will enable the features necessary to run WSL and install the Ubuntu distribution.
    :::

1. When prompted, restart Windows.
1. After restarting Windows, open the Ubuntu distribution and create a new Linux **User Name** and **Password** when prompted.

    :::{note}
    The **User Name** and **Password** is specific to each Linux distribution that you install and has no bearing on your Windows user name.
    :::

See [Set up a WSL development environment](https://learn.microsoft.com/en-us/windows/wsl/setup/environment) for more about installing WSL.
