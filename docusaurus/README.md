
# Docusaurus Development Notes (from owlang, 220707)

## 1. Set-up Node JS
Use a conda environment to set-up. If you don't have conda, check [this guide](https://medium.com/ayuth/install-anaconda-on-macos-with-homebrew-c94437d63a37).

```
$ conda update -n docusaurus -c conda-forge nodejs
$ conda activate docusaurus
```

## 2. Install Packages
If you don't have node packages already setup and installed in your docusaurus directory (`node_modules`), then you need to install them with node package manager:

```
cd scriptmanager/docusaurus
$ npm install
```
Sometimes an update is needed:

```
$ npm update
$ npm install
```

## 3. Develop
As you make changes to the source and markdown files, you can check your progress by running the site locally:
```
$ npm run start
```

## 4. Deploy
Before deploying the website, you need to check that the site will build properly. Test the website build by running:
```
$ npm run build
```
...and then test the website serving by running:
```
$ npm run serve
```

If the above tests run successfully, you can then deploy the website to update the `gh-pages` branch.

```
$ USE_SSH=true
$ GIT_USER=username npm run deploy
...
...
Password for 'https://username@github.com': <Your GitHub pass token>
...
...
Website is live at "https://CEGRcode.github.io/GenoPipe/".
```

# Website (initialization directions, old)

This website is built using [Docusaurus 2](https://docusaurus.io/), a modern static website generator.

### Installation

```
$ yarn
```

### Local Development

```
$ yarn start
```

This command starts a local development server and opens up a browser window. Most changes are reflected live without having to restart the server.

### Build

```
$ yarn build
```

This command generates static content into the `build` directory and can be served using any static contents hosting service.

### Deployment

Using SSH:

```
$ USE_SSH=true yarn deploy
```

Not using SSH:

```
$ GIT_USER=<Your GitHub username> yarn deploy
```

If you are using GitHub pages for hosting, this command is a convenient way to build the website and push to the `gh-pages` branch.
