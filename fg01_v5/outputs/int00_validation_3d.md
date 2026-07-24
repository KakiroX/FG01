**int00_validation_3d** — Closed-form |w| and retrograde efficiency checked against explicit 3-D inertial states at a common node.

| theta_deg | w_3d_m_s  | w_closed_form | rel_err_w | cos_psi_3d | cos_psi_closed_form | rel_err_cos_psi |
|-----------|-----------|---------------|-----------|------------|---------------------|-----------------|
| 0.1       | 12.92     | 12.92         | 0         | 0.0008727  | 0.0008727           | 2.77e-11        |
| 1         | 129.2     | 129.2         | 0         | 0.008727   | 0.008727            | 3.532e-13       |
| 8.5       | 1097      | 1097          | 2.072e-16 | 0.07411    | 0.07411             | 4.494e-15       |
| 45        | 5667      | 5667          | 1.605e-16 | 0.3827     | 0.3827              | 5.802e-16       |
| 90        | 1.047e+04 | 1.047e+04     | 0         | 0.7071     | 0.7071              | 0               |