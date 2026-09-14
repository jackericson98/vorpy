"""Temporary, serial call-boundary instrumentation for the CLI benchmark."""
from contextlib import contextmanager
from functools import wraps
import time


class Recorder:
    def __init__(self, wall=time.perf_counter, cpu=time.process_time):
        self.wall = wall
        self.cpu = cpu
        self.events = []
        self.stack = []

    @contextmanager
    def span(self, name, category):
        event = {
            'id': len(self.events), 'parent': self.stack[-1] if self.stack else None,
            'name': name, 'category': category, 'calls': 1,
            'wall_start': self.wall(), 'cpu_start': self.cpu(),
        }
        self.events.append(event)
        self.stack.append(event['id'])
        try:
            yield
        finally:
            event['cpu_end'] = self.cpu()
            event['wall_end'] = self.wall()
            event['wall_seconds'] = event['wall_end'] - event['wall_start']
            event['cpu_seconds'] = event['cpu_end'] - event['cpu_start']
            self.stack.pop()

    def report(self):
        events = [dict(event) for event in self.events]
        categories = {}
        for event in events:
            children = [child for child in events if child['parent'] == event['id']]
            category = categories.setdefault(event['category'], {
                'calls': 0, 'wall_seconds': 0.0, 'cpu_seconds': 0.0,
            })
            category['calls'] += 1
            for clock in ('wall', 'cpu'):
                child_time = sum(child[f'{clock}_seconds'] for child in children)
                remainder = event[f'{clock}_seconds'] - child_time
                if remainder < -1e-9:
                    raise ValueError(f'Negative {clock} remainder: {event["name"]}')
                event[f'{clock}_children_seconds'] = child_time
                event[f'{clock}_unattributed_seconds'] = remainder
                category[f'{clock}_seconds'] += remainder
        return {'events': events, 'exclusive_categories': categories}


@contextmanager
def instrument(recorder, hooks):
    """Restore every patched function, including when a timed call raises."""
    originals = []
    try:
        for owner, attribute, name, category in hooks:
            original = getattr(owner, attribute)

            def wrap(function, span_name, span_category):
                @wraps(function)
                def measured(*args, **kwargs):
                    with recorder.span(span_name, span_category):
                        return function(*args, **kwargs)
                return measured

            originals.append((owner, attribute, original))
            setattr(owner, attribute, wrap(original, name, category))
        yield
    finally:
        for owner, attribute, original in reversed(originals):
            setattr(owner, attribute, original)
