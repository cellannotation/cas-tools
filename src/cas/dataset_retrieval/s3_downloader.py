from cas.dataset_retrieval.dataset_retriever import DatasetRetriever


class S3Downloader(DatasetRetriever):
    def download_data(self):
        raise NotImplemented
