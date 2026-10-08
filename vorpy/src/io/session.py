"""Scientific session document; no GUI or calculation dependencies."""
from dataclasses import dataclass, field
from uuid import uuid4


@dataclass
class Session:
    systems: list = field(default_factory=list)
    networks: list = field(default_factory=list)
    workbench: dict = field(default_factory=dict)
    active_network_id: str | None = None
    archive_id: str = field(default_factory=lambda: str(uuid4()))

    @classmethod
    def from_network(cls, network):
        """Wrap only state that exists; never calculate missing analyses."""
        from .scientific_adapters import stable_id
        systems = [network.sys]
        networks = []
        for owner in (network.sys.groups or []) + (network.sys.ifaces or []):
            net = getattr(owner, 'net', None)
            if net is not None and all(net is not item for item in networks):
                networks.append(net)
        if all(network is not item for item in networks):
            networks.insert(0, network)
        return cls(systems, networks, active_network_id=stable_id(network))

    @property
    def active_network(self):
        return next((net for net in self.networks
                     if getattr(net, 'archive_id', None) == self.active_network_id),
                    self.networks[0] if self.networks else None)

    def capabilities(self, **filters):
        from .capabilities import inventory
        return inventory(self, **filters)
