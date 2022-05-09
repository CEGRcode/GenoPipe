import React from 'react';
import clsx from 'clsx';
import styles from './styles.module.css';

const FeatureList = [
  {
    title: 'EpitopeID',
    Svg: require('@site/static/genopipe-img/Figure1.svg').default,
    description: (
      <>
        Check your FASTQ files to confirm epitope sequence insertions in genomic DNA.
      </>
    ),
  },
  {
    title: 'DeletionID',
    Svg: require('@site/static/genopipe-img/Figure1.svg').default,
    description: (
      <>
        Check your BAM files for regions with significant depletions of reads.
      </>
    ),
  },
  {
    title: 'StrainID',
    Svg: require('@site/static/genopipe-img/Figure1.svg').default,
    description: (
      <>
        Check your BAM files for variants and see which strain your sample most closely resembles.
      </>
    ),
  },
];

function Feature({Svg, title, description}) {
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
