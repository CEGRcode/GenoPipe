// @ts-check
// Note: type annotations allow type checking and IDEs autocompletion

const lightCodeTheme = require('prism-react-renderer/themes/github');
const darkCodeTheme = require('prism-react-renderer/themes/dracula');

/** @type {import('@docusaurus/types').Config} */
const config = {
  title: 'GenoPipe',
  tagline: 'Run a background check on your NGS data',
  url: 'https://github.com/CEGRcode',
  baseUrl: '/GenoPipe/',
  onBrokenLinks: 'throw',
  onBrokenMarkdownLinks: 'warn',
  favicon: 'genopipe-img/logo.png',

  // GitHub pages deployment config.
  // If you aren't using GitHub pages, you don't need these.
  organizationName: 'CEGRcode', // Usually your GitHub org/user name.
  projectName: 'GenoPipe', // Usually your repo name.

  // Even if you don't use internalization, you can use this field to set useful
  // metadata like html lang. For example, if your site is Chinese, you may want
  // to replace "en" with "zh-Hans".
  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },

  presets: [
    [
      'classic',
      /** @type {import('@docusaurus/preset-classic').Options} */
      ({
        docs: {
          sidebarPath: require.resolve('./sidebars.js'),
          // Please change this to your repo.
          // Remove this to remove the "edit this page" links.
          editUrl:
            'https://github.com/CEGRcode/GenoPipe',
        },
        blog: {
          showReadingTime: true,
          // Please change this to your repo.
          // Remove this to remove the "edit this page" links.
          editUrl:
            'https://github.com/CEGRcode/GenoPipe',
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      }),
    ],
  ],

  themeConfig:
    /** @type {import('@docusaurus/preset-classic').ThemeConfig} */
    ({
      navbar: {
        title: 'GenoPipe',
        logo: {
          alt: 'GenoPipe Logo',
          src: 'genopipe-img/logo.png',
        },
        items: [
          {
            type: 'doc',
            docId: 'welcome',
            position: 'left',
            label: 'General',
          },
	        {
	          type: 'doc',
	          docId: 'epitopeid',
	          position: 'left',
	          label: 'EpitopeID',
	        },
	        {
	          type: 'doc',
	          docId: 'deletionid',
	          position: 'left',
	          label: 'DeletionID',
	        },
	        {
	          type: 'doc',
	          docId: 'strainid',
	          position: 'left',
	          label: 'StrainID',
	        },
          {
            href: 'https://github.com/CEGRcode/GenoPipe',
            label: 'GitHub',
            position: 'right',
          },
        ],
      },
      footer: {
        style: 'dark',
        links: [
          {
            title: 'Module Docs',
            items: [
              {
                label: 'EpitopeID',
                to: '/docs/epitopeid',
              },
              {
                label: 'DeletionID',
                to: '/docs/deletionid',
              },
              {
                label: 'StrainID',
                to: '/docs/strainid',
              },
            ],
          },
          {
            title: 'Community',
            items: [
              {
                label: 'Pugh Lab Website',
                href: 'https://pughlab.mbg.cornell.edu/',
              },
            ],
          },
          {
            title: 'More',
            items: [
              {
                label: 'GitHub',
                href: 'https://github.com/CEGRcode/GenoPipe',
              },
            ],
          },
        ],
        copyright: `Copyright © ${new Date().getFullYear()} GenoPipe. Built with Docusaurus.`,
      },
      prism: {
        theme: lightCodeTheme,
        darkTheme: darkCodeTheme,
      },
    }),
};

module.exports = config;
