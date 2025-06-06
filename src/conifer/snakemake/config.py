from nbis.config import Schema
from ruamel.yaml import YAML

from conifer.cli import PKG_DIR

WORKFLOW_ROOT = PKG_DIR / "workflow"
SNAKEMAKE_ROOT = WORKFLOW_ROOT / "snakemake"


class SchemaFiles:
    CONFIGURATION_SCHEMA = WORKFLOW_ROOT / "schemas" / "config.schema.yaml"
    SAMPLES_SCHEMA = WORKFLOW_ROOT / "schemas" / "samples.schema.yaml"


def get_schema(schema="CONFIGURATION_SCHEMA"):
    schemafile = getattr(SchemaFiles, schema)
    with open(schemafile) as fh:
        schema = YAML().load(fh)
    return Schema(schema)
