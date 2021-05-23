class Diana: public Algorithm{
    protected:
        void setupArguments();
        void performClustering() const;
        void fetchResult() const;
        virtual void create_dm() const;
        virtual void initialization() const;
        virtual void division() const;
        virtual void do_split(Size ind) const;
        virtual void create_cluster(const std::set<Size> ele, Size ind) const;
        mutable iirMapA _dm;
        mutable std::set<Size> _unsplitClusters;
        mutable std::map<Size, Real> _clusterDiameter;
        mutable std::map<Size, boost::shared_ptr<LeafNode>> _leaf;
        mutable std::map<Size, boost::shared_ptr<InternalNode>> _internal;
        mutable std::set<Size> _clusterID;
        boost::shared_ptr<Distance> _distance;
}

void Diana::initialization() const {
    Size n = _ds->size();
    Size id = 2*n - 2;
    boost::shared_ptr<InternalNode> pin(new InternalNode(id));
    
}

