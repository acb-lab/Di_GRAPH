# Config models

Pydantic v2 models for Di-GRAPH configuration. All models are validated on startup via [`load_config()`][digraph.config.load_config].

---

::: digraph.config
    options:
      members:
        - load_config
        - DiGraphConfig
        - PathsConfig
        - GenomeConfig
        - MATCoordinates
        - PolymorphismConfig
        - TrimConfig
        - ExperimentConfig
        - SampleConfig
        - ResourceConfig
      show_source: false
      show_root_heading: true
