from contextlib import contextmanager
import email
import re
from time import sleep
from urllib.error import HTTPError, URLError
from collections.abc import Generator, Iterator
from typing import IO

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature


def get_accessions_from_search(taxid: str, min_len: int, max_len: int, email: str, api_key: str | None = None) -> list[str]:
    term = f"txid{taxid}[Organism] AND {min_len}:{max_len}[SLEN]"

    # use the credentials only within this context to avoid affecting global state
    with temporary_entrez_credentials(email=email, api_key=api_key):

        with Entrez.esearch(db="nuccore", term=term, usehistory="y") as handle:
            search_results = Entrez.read(handle)
        
        web_env = search_results["WebEnv"]
        query_key = search_results["QueryKey"]
        count = int(search_results["Count"])
        
        accessions = []
        # NCBI caps esummary at 100 results per request, so paginate in chunks
        for start in range(0, count, 100):
            with Entrez.esummary(
                db="nuccore",
                WebEnv=web_env,
                query_key=query_key,
                retstart=start,
                retmax=100
            ) as handle:
                summaries = Entrez.read(handle)
                accessions.extend([doc["AccessionVersion"] for doc in summaries])
                
        return accessions


def extract_feature(
    record: SeqRecord,
    feature_type: str,
    qualifier: str = "gene",
    value: str | None = None,
) -> SeqFeature | None:
    """
    Return the first feature matching the requested criteria.

    Parameters
    ----------
    record
        GenBank SeqRecord.
    feature_type
        Feature type (e.g. "CDS", "gene", "tRNA").
    qualifier
        Qualifier used for matching (default: "gene").
    value
        Desired qualifier value. If None, return the first feature of
        the requested type.

    Returns
    -------
    SeqFeature or None
    """
    for feature in record.features:
        if feature.type != feature_type:
            continue

        if value is None:
            return feature

        if value in feature.qualifiers.get(qualifier, []):
            return feature

    return None


@contextmanager
def entrez_efetch_with_retries(
    *,
    db: str,
    ids: list[str],
    rettype: str,
    retmode: str,
    retries: int,
) -> Generator[IO[str], None, None]:
    """Yield an Entrez efetch handle with exponential-backoff retries."""
    for attempt in range(retries):
        try:
            with Entrez.efetch(
                db=db,
                id=",".join(ids),
                rettype=rettype,
                retmode=retmode,
            ) as handle:
                yield handle
            return
        except (HTTPError, URLError):
            if attempt == retries - 1:
                raise
            sleep(2 ** attempt)

@contextmanager
def temporary_entrez_credentials(email: str, api_key: str | None = None):
    original_email = Entrez.email
    original_api_key = Entrez.api_key
    
    try:
        if email:
            Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
        yield
    finally:
        # Restore original values
        Entrez.email = original_email
        Entrez.api_key = original_api_key


def fetch_genbank_records(
    accessions: list[str],
    email: str,
    api_key: str | None = None,
    batch_size: int = 100,
    retries: int = 3,
    feature_type: str | None = None,
    qualifier: str = "gene",
    value: str | None = None,
) -> Iterator[SeqRecord]:
    """
    Fetch GenBank SeqRecords in batches.

    Parameters
    ----------
    accessions
        NCBI nucleotide accession identifiers to fetch.
    email
        Contact email required by NCBI Entrez.
    api_key
        Optional NCBI API key.
    batch_size
        Number of accessions to fetch per request.
    retries
        Number of retries for transient HTTP/network errors.
    feature_type
        Optional feature type used to filter records (e.g. "CDS", "gene").
    qualifier
        Feature qualifier used when filtering by ``value``.
    value
        Optional qualifier value used for feature filtering.

    Yields
    ------
    SeqRecord
        If ``feature_type`` is not provided, yields full fetched records.
        If ``feature_type`` is provided, yields new SeqRecord objects with
        sequences extracted from the matched feature.

    Notes
    -----
    This generator never yields SeqFeature objects.
    """

    with temporary_entrez_credentials(email=email, api_key=api_key):

        def build_record_id(record_id: str, ft: str) -> str:
            """Build a deterministic output record ID from active feature filters."""
            parts = (("ft", ft), ("q", qualifier), ("v", value))
            suffix = "_".join(
                re.sub(r"[^A-Za-z0-9._-]+", "-", f"{key}-{part}").strip("-")
                for key, part in parts
                if part
            )
            return f"{record_id}_{suffix}"

        for start in range(0, len(accessions), batch_size):
            batch = accessions[start:start + batch_size]

            with entrez_efetch_with_retries(
                db="nucleotide",
                ids=batch,
                rettype="gb",
                retmode="text",
                retries=retries,
            ) as handle:

                for record in SeqIO.parse(handle, "genbank"):

                    if feature_type is None:
                        yield record
                    else:
                        feature = extract_feature(
                            record,
                            feature_type=feature_type,
                            qualifier=qualifier,
                            value=value,
                        )
                        if feature is not None:
                            yield SeqRecord(
                                seq=feature.extract(record.seq),
                                id=build_record_id(record.id, feature_type),
                                name=record.name,
                                description=record.description,
                            )