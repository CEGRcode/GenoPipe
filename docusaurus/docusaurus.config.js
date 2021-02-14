module.exports = {
  title: 'GenoPipe',
  tagline: 'Run a background check on your NGS data',
  url: 'https://your-docusaurus-test-site.com',
  baseUrl: '/',
  onBrokenLinks: 'throw',
  onBrokenMarkdownLinks: 'warn',
  favicon: 'genopipe-img/logo.png',
  organizationName: 'facebook', // Usually your GitHub org/user name.
  projectName: 'docusaurus', // Usually your repo name.
  themeConfig: {
    navbar: {
      title: 'GenoPipe',
      logo: {
        alt: 'GenoPipe Logo',
        src: 'genopipe-img/logo.png',
      },
      items: [
        {
          to: 'docs/welcome',
          activeBasePath: 'docs',
          label: 'General',
          position: 'left',
        },
        {
          to: 'docs/epitopeid',
          label: 'EpitopeID',
          position: 'left',
        },
        {
          to: 'docs/deletionid',
          label: 'DeletionID',
          position: 'left',
        },
        {
          to: 'docs/strainid',
          label: 'StrainID',
          position: 'left',
        },
        {to: 'blog', label: 'Blog', position: 'left'},
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
              to: 'docs/epitopeid',
            },
            {
              label: 'DeletionID',
              to: 'docs/deletionid',
            },
            {
              label: 'StrainID',
              to: 'docs/strainid',
            },
          ],
        },
        {
          title: 'Community',
          items: [
            {
              label: 'Blog',
              to: 'blog',
            },
            {
              label: 'GitHub',
              href: 'https://github.com/CEGRcode/GenoPipe',
            },
          ],
        },
      ],
      copyright: `Copyright © ${new Date().getFullYear()} My Project, Inc. Built with Docusaurus.`,
    },
  },
  presets: [
    [
      '@docusaurus/preset-classic',
      {
        docs: {
          sidebarPath: require.resolve('./sidebars.js'),
          // Please change this to your repo.
          editUrl:
            'https://github.com/CEGRcode/GenoPipe',
        },
        blog: {
          showReadingTime: true,
          // Please change this to your repo.
          editUrl:
            'https://github.com/CEGRcode/GenoPipe',
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      },
    ],
  ],
};
