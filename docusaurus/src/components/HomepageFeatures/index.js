import React from "react";
import clsx from "clsx";
import styles from "./styles.module.css";
import { translate } from "@docusaurus/Translate";

const FeatureList = [
  {
    title: "EpitopeID",
    Svg: require("@site/static/genopipe-img/EpitopeID_icon.svg").default,
    description: (
      <>
        Check your FASTQ files to confirm epitope sequence insertions in genomic
        DNA.
      </>
    ),
    link: "docs/epitopeid",
  },
  {
    title: "DeletionID",
    Svg: require("@site/static/genopipe-img/DeletionID_icon.svg").default,
    description: (
      <>
        Check your BAM files for regions with significant depletions of reads.
      </>
    ),
    link: "docs/deletionid",
  },
  {
    title: "StrainID",
    Svg: require("@site/static/genopipe-img/StrainID_icon.svg").default,
    description: (
      <>
        Check your BAM files for variants and see which strain your sample most
        closely resembles.
      </>
    ),
    link: "docs/strainid",
  },
];

/*function Feature({Svg, title, description}) {
  return (
    <div className={clsx('col col--4')}>
      <div className="text--center">
        <Svg className={styles.featureSvg} role="img" />
      </div>
      <div className="text--center padding-horiz--md">
        <h3>{title}</h3>
        <p>{description}</p>
      </div>
    </div>
  );
}*/

function Feature({ Svg, title, description, link }) {
  return (
    <div
      class="card-demo col col--4 text--center"
      style={{ paddingTop: 20 + "px", paddingBottom: 80 + "px" }}
    >
      <div
        class="card"
        style={{
          paddingTop: 25 + "px",
          paddingBottom: 25 + "px",
          boxShadow: "10px 10px 30px lightgrey",
        }}
      >
        <div class="card__image">
          <Svg className={styles.featureSvg} role="img" />
        </div>
        <div class="card__body">
          <h3>{title}</h3>
          <p>{description}</p>
        </div>
        <div class="card__footer">
          <a class="" href={link}>
            Read More
            <svg
              width="13.5"
              height="13.5"
              aria-hidden="true"
              viewBox="0 0 24 24"
              class="iconExternalLink_node_modules-@docusaurus-theme-classic-lib-next-theme-IconExternalLink-styles-module"
            >
              <path
                fill="currentColor"
                d="M21 13v10h-21v-19h12v2h-10v15h17v-8h2zm3-12h-10.988l4.035 4-6.977 7.07 2.828 2.828 6.977-7.07 4.125 4.172v-11z"
              ></path>
            </svg>
          </a>
        </div>
      </div>
    </div>
  );
}

export default function HomepageFeatures() {
  return (
    <section className={styles.features}>
      <div className="container">
        <div className="row">
          {FeatureList.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}
