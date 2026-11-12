{{ fullname | escape | underline }}

{% if modules %}
.. automodule:: {{ fullname }}
   :no-members:

.. rubric:: Modules

.. autosummary::
   :toctree:
   :recursive:
{% for item in modules %}
   {{ item }}
{%- endfor %}
{% else %}
.. automodule:: {{ fullname }}
   :members:
{% endif %}
